// adapted from schreier.c of nauty 2.5R9

/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include "schreier_sequential.h"

static permnode id_permnode;
/* represents identity, no actual content, doesn't need TLS_ATTR */
#define ID_PERMNODE (&id_permnode)

#if !MAXN
DYNALLSTAT(int,workperm,workperm_sz);
DYNALLSTAT(int,workperm2,workperm2_sz);
DYNALLSTAT(int,workpermA,workpermA_sz);
DYNALLSTAT(int,workpermB,workpermB_sz);
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(set,workset2,workset2_sz);
#else
static TLS_ATTR int workperm[MAXN];
static TLS_ATTR int workperm2[MAXN];
static TLS_ATTR int workpermA[MAXN];
static TLS_ATTR int workpermB[MAXN];
static TLS_ATTR set workset[MAXM];
static TLS_ATTR set workset2[MAXM];
#endif

static TLS_ATTR _schreier *_schreier_freelist = NULL;
/* Freelist of scheier structures connected by next field.
     * vec, pwr and orbits fields are assumed allocated. */
static TLS_ATTR permnode *_permnode_freelist = NULL;
/* Freelist of permnode structures connected by next field.
     * p[] is assumed extended. */

static TLS_ATTR int schreierfails = SCHREIERFAILS;

#define TMP

static boolean _filterschreier(_schreier*,int*,permnode**,boolean,int,int,int*);
#define PNCODE(x) ((int)(((size_t)(x)>>3)&0xFFFUL))

/* #define TESTP(id,p,n) testispermutation(id,p,n) */
#define TESTP(id,p,n)

/************************************************************************/

static void
_testispermutation(int id, int *p, int n)
/* For debugging purposes, crash with a message if p[0..n-1] is
   not a permutation. */
{
    int i,m;
    DYNALLSTAT(set,seen,seen_sz);

    for (i = 0; i < n; ++i)
        if (p[i] < 0 || p[i] > n) break;

    if (i < n)
    {
        fprintf(stderr,">E Bad permutation (id=%d): n=%d p[%d]=%d\n",
                id,n,i,p[i]);
        exit(1);
    }

    m = SETWORDSNEEDED(n);
    DYNALLOC1(set,seen,seen_sz,m,"malloc seen");
    EMPTYSET(seen,m);

    for (i = 0; i < n; ++i)
    {
        if (ISELEMENT(seen,p[i]))
        {
            fprintf(stderr,
                    ">E Bad permutation (id=%d): n=%d p[%d]=%d is a repeat\n",
                    id,n,i,p[i]);
            exit(1);
        }
                ADDELEMENT(seen,p[i]);
    }
}

/************************************************************************/

int
_schreier_fails(int nfails)
/* Set the number of consecutive failures for filtering;
 * A value of <= 0 defaults to SCHREIERFAILS.
 * The function value is the previous setting. */
{
    int prev;

    prev = schreierfails;

    if (nfails <= 0) schreierfails = SCHREIERFAILS;
    else             schreierfails = nfails;

    return prev;
}

/************************************************************************/

static void
_clearfreelists(void)
/* Clear the schreier and permnode freelists */
{
    _schreier *sh,*nextsh;
    permnode *p,*nextp;

    nextsh = _schreier_freelist;
    while (nextsh)
    {
        sh = nextsh;
        nextsh = sh->next;
        free(sh->vec);
        free(sh->pwr);
        free(sh->orbits);
        free(sh->orbits_sz);
        free(sh);
    }
    _schreier_freelist = NULL;

    nextp = _permnode_freelist;
    while (nextp)
    {
        p = nextp;
        nextp = p->next;
        free(p);
    }
    _permnode_freelist = NULL;
}

/************************************************************************/

static permnode
*_newpermnode(int n)
/* Allocate a new permode structure, with initialized next/prev fields */
{
    permnode *p;

    while (_permnode_freelist)
    {
        p = _permnode_freelist;
        _permnode_freelist = p->next;
        if (p->nalloc >= n && p->nalloc <= n+100)
        {
            p->next = p->prev = NULL;
            p->mark = 0;
            return p;
        }
        else
            free(p);
    }

    p = (permnode*) malloc(sizeof(permnode)+(n-2)*sizeof(int));

    if (p == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newpermnode()\n");
        exit(1);
    }

    p->next = p->prev = NULL;
    p->nalloc = n;

    return p;
}

/************************************************************************/

static _schreier
*_newschreier(int n)
/* Allocate a new schreier structure, with initialised next field */
{
    _schreier *sh;

    while (_schreier_freelist)
    {
        sh = _schreier_freelist;
        _schreier_freelist = sh->next;
        if (sh->nalloc >= n && sh->nalloc <= n+100)
        {
            sh->next = NULL;
            return sh;
        }
        else
        {
            free(sh->vec);
            free(sh->pwr);
            free(sh->orbits);
            free(sh->orbits_sz);
            free(sh);
        }
    }

    sh = (_schreier*) malloc(sizeof(_schreier));

    if (sh == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->vec = (permnode**) malloc(sizeof(permnode*)*n);
    sh->pwr = (int*) malloc(sizeof(int)*n);
    sh->orbits = (int*) malloc(sizeof(int)*n);
    sh->orbits_sz = (int*) malloc(sizeof(int)*n);

    if (sh->vec == NULL || sh->pwr == NULL || sh->orbits == NULL || sh->orbits_sz == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->next = NULL;
    sh->nalloc = n;

    return sh;
}

/************************************************************************/

void
_freeschreier(_schreier **gp, permnode **gens)
/* Free schreier structure and permutation ring.  Assume this is everything. */
/* Use NULL for arguments which don't need freeing. */
{
    _schreier *sh;
    _schreier *nextsh;
    permnode *p,*nextp;

    if (gp && *gp)
    {
        nextsh = *gp;
        while (nextsh)
        {
            sh = nextsh;
            nextsh = sh->next;
            sh->next = _schreier_freelist;
            _schreier_freelist = sh;
        }
        *gp = NULL;
    }

    if (gens && *gens)
    {
        p = *gens;
        do
        {
            nextp = p->next;
            p->next = _permnode_freelist;
            _permnode_freelist = p;
            p = nextp;
        } while (p != *gens);
        *gens = NULL;
    }
}

/************************************************************************/

permnode*
_findpermutation(permnode *pn, int *p, int n)
/* Return a pointer to permutation p in the circular list,
 * or NULL if it isn't present. */
{
    permnode *rn;
    int i;

    if (!pn) return NULL;

    rn = pn;
    do
    {
        for (i = 0; i < n; ++i)
            if (rn->p[i] != p[i]) break;
        if (i == n) return rn;
        rn = rn->next;
    } while (rn != pn);

    return NULL;
}

/************************************************************************/

void
_addpermutation(permnode **ring, int *p, int n)
/* Add new permutation to circular list, marked.
 * and return pointer to it in *ring. */
{
    permnode *pn,*rn;

    pn = _newpermnode(n);
    rn = *ring;

    memcpy(pn->p,p,n*sizeof(int));

    if (!rn)
        pn->next = pn->prev = pn;
    else
    {
        pn->next = rn->next;
        pn->prev = rn;
        rn->next = pn->next->prev = pn;
    }

    pn->refcount = 0;
    pn->mark = 1;
    *ring = pn;
}

/************************************************************************/

static void
_addpermutationunmarked(permnode **ring, int *p, int n)
/* Add new permutation to circular list, not marked.
 * and return pointer to it in *ring. */
{
    TESTP(3,p,n);
    _addpermutation(ring,p,n);
    (*ring)->mark = 0;
}

/************************************************************************/

boolean
_addgenerator(_schreier **gp, permnode **ring, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    TESTP(2,p,n);
    return _filterschreier(*gp,p,ring,FALSE,-1,n,NULL);
}

/************************************************************************/

static void
_delpermnode(permnode **ring)
/* Delete permnode at head of circular list, making the next node head. */
{
    permnode *newring;

    if (!*ring) return;

    if ((*ring)->next == *ring)
        newring = NULL;
    else
    {
        newring = (*ring)->next;
        newring->prev = (*ring)->prev;
        (*ring)->prev->next = newring;
    }

    (*ring)->next = _permnode_freelist;
    _permnode_freelist = *ring;

    *ring = newring;
}

/************************************************************************/

void
_deleteunmarked(permnode **ring)
/* Delete all permutations in the ring that are not marked */
{
    permnode *pn,*firstmarked;

    pn = *ring;
    firstmarked = NULL;

    while (pn != NULL && pn != firstmarked)
    {
        if (pn->mark)
        {
            if (!firstmarked) firstmarked = pn;
            pn = pn->next;
        }
        else
            _delpermnode(&pn);
    }

    *ring = pn;
}

/************************************************************************/

static void
_clearvector(permnode **vec, permnode **ring, int n)
/* clear vec[0..n-1], freeing permnodes that have no other references
 * and are not marked */
{
    int i;

    for (i = 0; i < n; ++i)
        if (vec[i])
        {
            if (vec[i] != ID_PERMNODE)
            {
                --(vec[i]->refcount);
                if (vec[i]->refcount == 0 && !vec[i]->mark)
                {
                    *ring = vec[i];
                    _delpermnode(ring);
                }
            }
            vec[i] = NULL;
        }
}

/************************************************************************/

static void
_initschreier(_schreier *sh, int n)
/* Initialise schreier structure to trivial orbits and empty vector */
{
    int i;

    sh->fixed = -1;
    for (i = 0; i < n; ++i)
    {
        sh->vec[i] = NULL;
        sh->orbits[i]    = i;
        sh->orbits_sz[i] = 1;
    }
}

/************************************************************************/

void
_newgroup(_schreier **sh, permnode **ring, int n)
/* Make the trivial group, allow for ring to be set elsewhere */
{
    *sh = _newschreier(n);
    _initschreier(*sh,n);
    if (ring) *ring = NULL;
}

/************************************************************************/

static void
_applyperm(int *wp, int *p, int k, int n)
/* Apply the permutation p, k times to each element of wp */
{
    int i,j,cyclen,kk,m;

    TESTP(1,p,n);

    if (k <= 5)
    {
        if (k == 0)
            return;
        else if (k == 1)
            for (i = 0; i < n; ++i) wp[i] = p[wp[i]];
        else if (k == 2)
            for (i = 0; i < n; ++i) wp[i] = p[p[wp[i]]];
        else if (k == 3)
            for (i = 0; i < n; ++i) wp[i] = p[p[p[wp[i]]]];
        else if (k == 4)
            for (i = 0; i < n; ++i) wp[i] = p[p[p[p[wp[i]]]]];
        else if (k == 5)
            for (i = 0; i < n; ++i) wp[i] = p[p[p[p[p[wp[i]]]]]];
    }
    else if (k <= 19)
    {
#if !MAXN
        DYNALLOC1(int,workpermA,workpermA_sz,n,"applyperm");
#endif
        for (i = 0; i < n; ++i) workpermA[i] = p[p[p[i]]];
        for (; k >= 6; k -= 6)
            for (i = 0; i < n; ++i) wp[i] = workpermA[workpermA[wp[i]]];
        if (k == 1)
            for (i = 0; i < n; ++i) wp[i] = p[wp[i]];
        else if (k == 2)
            for (i = 0; i < n; ++i) wp[i] = p[p[wp[i]]];
        else if (k == 3)
            for (i = 0; i < n; ++i) wp[i] = workpermA[wp[i]];
        else if (k == 4)
            for (i = 0; i < n; ++i) wp[i] = p[workpermA[wp[i]]];
        else if (k == 5)
            for (i = 0; i < n; ++i) wp[i] = p[p[workpermA[wp[i]]]];
    }
    else
    {
        m = SETWORDSNEEDED(n);
#if !MAXN
        DYNALLOC1(int,workpermA,workpermA_sz,n,"applyperm");
        DYNALLOC1(int,workpermB,workpermB_sz,n,"applyperm");
        DYNALLOC1(set,workset2,workset2_sz,m,"applyperm");
#endif

        EMPTYSET(workset2,m);

        /* We will construct p^k in workpermB one cycle at a time. */

        for (i = 0; i < n; ++i)
        {
            if (ISELEMENT(workset2,i)) continue;
            if (p[i] == i)
                workpermB[i] = i;
            else
            {
                cyclen = 1;
                workpermA[0] = i;
                for (j = p[i]; j != i; j = p[j])
                {
                    workpermA[cyclen++] = j;
                            ADDELEMENT(workset2,j);
                }
                kk = k % cyclen;
                for (j = 0; j < cyclen; ++j)
                {
                    workpermB[workpermA[j]] = workpermA[kk];
                    if (++kk == cyclen) kk = 0;
                }
            }
        }
        for (i = 0; i < n; ++i) wp[i] = workpermB[wp[i]];
    }
}

/************************************************************************/

static boolean
_filterschreier(_schreier *gp, int *p, permnode **ring,
               boolean ingroup, int maxlevel, int n, int* minimal_base)
/* Filter permutation p up to level maxlevel of gp.
 * Use ingroup=TRUE if p is known to be in the group, otherwise
 * at least one equivalent generator is added unless it is proved
 * (nondeterministically) that it is in the group already. 
 * maxlevel < 0 means no limit, maxlevel=0 means top level only, etc.
 * Return TRUE iff some change is made. */
{
    int i,j,j1,j2,lev;
    int ipwr;
    _schreier *sh;
    int *orbits,*pwr;
    permnode **vec,*curr;
    boolean changed,lchanged,ident;
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"filterschreier");
#endif

    memcpy(workperm,p,n*sizeof(int));

    if (*ring && p == (*ring)->p)
    {
        ingroup = TRUE;
        curr = *ring;
    }
    else
        curr = NULL;

    /* curr is the location of workperm in ring, if anywhere */

    sh = gp;
    changed = FALSE;
    if (maxlevel < 0) maxlevel = n+1;

    for (lev = 0; lev <= maxlevel; ++lev)
    {
        for (i = 0; i < n; ++i) if (workperm[i] != i) break;
        ident = (i == n);
        if (ident) break;

        lchanged = FALSE;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;
        for (i = 0; i < n; ++i)
        {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[workperm[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];

            if (j1 != j2)
            {
                lchanged = TRUE;
                // prefer minimal base
                if ((j1 < j2 || j1 == minimal_base[lev]) && (j2 != minimal_base[lev])) {
                    orbits[j2] = j1;
                    sh->orbits_sz[j1] +=  sh->orbits_sz[j2];
                } else {
                    orbits[j1] = j2;
                    sh->orbits_sz[j2] += sh->orbits_sz[j1];
                }
            }
        }
        if (lchanged)
            for (i = 0; i < n; ++i) {
                if(orbits[i] != orbits[orbits[i]]) {
                    //sh->orbits_sz[orbits[orbits[i]]] += sh->orbits_sz[orbits[i]];
                    orbits[i] = orbits[orbits[i]];
                }
            }

        if (lchanged) changed = TRUE;

        if (sh->fixed >= 0)
        {
            for (i = 0; i < n; ++i)
                if (vec[i] && !vec[workperm[i]])
                {
                    changed = TRUE;
                    ipwr = 0;
                    for (j = workperm[i]; !vec[j] ; j = workperm[j]) ++ipwr;

                    for (j = workperm[i]; !vec[j] ; j = workperm[j])
                    {
                        if (!curr)
                        {
                            if (!ingroup) _addpermutation(ring,workperm,n);
                            else  _addpermutationunmarked(ring,workperm,n);
                            ingroup = TRUE;
                            curr = *ring;
                        }
                        vec[j] = curr;
                        pwr[j] = ipwr--;
                        ++curr->refcount;
                    }
                }

            j = workperm[sh->fixed];

            while (j != sh->fixed)
            {
                _applyperm(workperm,vec[j]->p,pwr[j],n);
                curr = NULL;
                j = workperm[sh->fixed];
            }
            sh = sh->next;
        }
        else
            break;
    }

    if (!ident && !ingroup)
    {
        changed = TRUE;
        _addpermutation(ring,p,n);
    }

    return changed;
}

/************************************************************************/

boolean
_expandschreier(_schreier *gp, permnode **ring, int n, int* minimal_base)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i,j,nfails,wordlen,skips;
    boolean changed;
    permnode *pn;
#if !MAXN
    DYNALLOC1(int,workperm2,workperm2_sz,n,"expandschreier");
#endif

    pn = *ring;
    if (pn == NULL) return FALSE;

    nfails = 0;
    changed = FALSE;

    for (skips = KRAN(17); --skips >= 0; ) pn = pn->next;

    memcpy(workperm2,pn->p,n*sizeof(int));

    while (nfails < schreierfails)
    {
        wordlen = 1 + KRAN(3);
        for (j = 0; j < wordlen; ++j)
        {
            for (skips = KRAN(17); --skips >= 0; ) pn = pn->next;
            for (i = 0; i < n; ++i) workperm2[i] = pn->p[workperm2[i]];
        }
        if (_filterschreier(gp,workperm2,ring,TRUE,-1,n,minimal_base))
        {
            changed = TRUE;
            nfails = 0;
        }
        else
            ++nfails;
    }

    return changed;
}

/************************************************************************/

int*
_getorbits(int *fix, int nfix, _schreier *gp, permnode **ring, int n, int* minimal_base, int** orbits_sz)
/* Get a pointer to the orbits for this partial base. The pointer
 * remains valid until pruneset(), getorbits(), getorbitsmin()
 * or grouporder() is called with an incompatible base (neither a
 * prefix nor an extension). The contents of the array pointed to
 * MUST NOT BE MODIFIED by the calling program.
 */
{
    int k;
    _schreier *sh,*sha;

    sh = gp;
    for (k = 0; k < nfix; ++k)
    {
        if (sh->fixed != fix[k]) break;
        sh = sh->next;
    }

    if (k == nfix) {
        *orbits_sz = sh->orbits_sz;
        return sh->orbits;
    }

    sh->fixed = fix[k];
    _clearvector(sh->vec,ring,n);
    sh->vec[fix[k]] = ID_PERMNODE;

    for (sha = sh->next; sha ; sha = sha->next) _clearvector(sha->vec,ring,n);

    for (++k; k <= nfix; ++k)
    {
        if (!sh->next) sh->next = _newschreier(n);
        sh = sh->next;
        _initschreier(sh,n);
        if (k < nfix)
        {
            sh->fixed = fix[k];
            sh->vec[fix[k]] = ID_PERMNODE;
        }
        else
            sh->fixed = -1;
    }

    if (*ring) _expandschreier(gp,ring,n,minimal_base);
    *orbits_sz = sh->orbits_sz;
    return sh->orbits;
}

/************************************************************************/

void
_schreier_freedyn(void)
{
#if !MAXN
    DYNFREE(workperm,workperm_sz);
    DYNFREE(workperm2,workperm2_sz);
    DYNFREE(workpermA,workpermA_sz);
    DYNFREE(workpermB,workpermB_sz);
    DYNFREE(workset,workset_sz);
    DYNFREE(workset2,workset2_sz);
#endif
    _clearfreelists();
}
