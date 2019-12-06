// adapted from schreier.c of nauty 2.5R9

/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include "schreier_sequential.h"

static sequential_permnode id_permnode;
#define ID_PERMNODE (&id_permnode)
DYNALLSTAT(int,workperm,workperm_sz);
DYNALLSTAT(int,workperm2,workperm2_sz);
DYNALLSTAT(int,workpermA,workpermA_sz);
DYNALLSTAT(int,workpermB,workpermB_sz);
DYNALLSTAT(set,workset,workset_sz);
DYNALLSTAT(set,workset2,workset2_sz);


static sequential_schreier *_schreier_freelist = NULL;
static sequential_permnode *_permnode_freelist = NULL;
static int schreierfails = SCHREIERFAILS;

static bool _filterschreier(sequential_schreier*, int*, sequential_permnode**, bool, int, int, int*);

/************************************************************************/

int _schreier_fails(int nfails) {
    int prev;

    prev = schreierfails;

    if (nfails <= 0) schreierfails = SCHREIERFAILS;
    else             schreierfails = nfails;

    return prev;
}

/************************************************************************/

static void _clearfreelists(void) {
    sequential_schreier *sh,*nextsh;
    sequential_permnode *p,*nextp;

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

static sequential_permnode *_newpermnode(int n) {
    sequential_permnode *p;

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

    p = (sequential_permnode*) malloc(sizeof(sequential_permnode) + (n - 2) * sizeof(int));

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

static sequential_schreier *_newschreier(int n) {
    sequential_schreier *sh;

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

    sh = (sequential_schreier*) malloc(sizeof(sequential_schreier));

    if (sh == NULL)
    {
        fprintf(ERRFILE,">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->vec = (sequential_permnode**) malloc(sizeof(sequential_permnode*) * n);
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

void _freeschreier(sequential_schreier **gp, sequential_permnode **gens) {
    sequential_schreier *sh;
    sequential_schreier *nextsh;
    sequential_permnode *p,*nextp;

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

sequential_permnode* _findpermutation(sequential_permnode *pn, int *p, int n) {
    sequential_permnode *rn;
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

void _addpermutation(sequential_permnode **ring, int *p, int n) {
    sequential_permnode *pn,*rn;

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

static void _addpermutationunmarked(sequential_permnode **ring, int *p, int n) {
    _addpermutation(ring,p,n);
    (*ring)->mark = 0;
}

/************************************************************************/

bool _addgenerator(sequential_schreier **gp, sequential_permnode **ring, int *p, int n) {
    return _filterschreier(*gp,p,ring,false,-1,n,NULL);
}

/************************************************************************/

static void _delpermnode(sequential_permnode **ring) {
    sequential_permnode *newring;

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

void _deleteunmarked(sequential_permnode **ring) {
    sequential_permnode *pn,*firstmarked;

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

static void _clearvector(sequential_permnode **vec, sequential_permnode **ring, int n) {
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

static void _initschreier(sequential_schreier *sh, int n) {
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

void _newgroup(sequential_schreier **sh, sequential_permnode **ring, int n) {
    *sh = _newschreier(n);
    _initschreier(*sh,n);
    if (ring) *ring = NULL;
}

/************************************************************************/

static void _applyperm(int *wp, int *p, int k, int n) {
    int i,j,cyclen,kk,m;

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

static bool _filterschreier(sequential_schreier *gp, int *p, sequential_permnode **ring,
                            bool ingroup, int maxlevel, int n, int* minimal_base) {
    int i,j,j1,j2,lev;
    int ipwr;
    sequential_schreier *sh;
    int *orbits,*pwr;
    sequential_permnode **vec,*curr;
    bool changed,lchanged,ident;
#if !MAXN
    DYNALLOC1(int,workperm,workperm_sz,n,"filterschreier");
#endif

    memcpy(workperm,p,n*sizeof(int));

    if (*ring && p == (*ring)->p)
    {
        ingroup = true;
        curr = *ring;
    }
    else
        curr = NULL;

    /* curr is the location of workperm in ring, if anywhere */

    sh = gp;
    changed = false;
    if (maxlevel < 0) maxlevel = n+1;

    for (lev = 0; lev <= maxlevel; ++lev)
    {
        for (i = 0; i < n; ++i) if (workperm[i] != i) break;
        ident = (i == n);
        if (ident) break;

        lchanged = false;
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
                lchanged = true;
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

        if (lchanged) changed = true;

        if (sh->fixed >= 0)
        {
            for (i = 0; i < n; ++i)
                if (vec[i] && !vec[workperm[i]])
                {
                    changed = true;
                    ipwr = 0;
                    for (j = workperm[i]; !vec[j] ; j = workperm[j]) ++ipwr;

                    for (j = workperm[i]; !vec[j] ; j = workperm[j])
                    {
                        if (!curr)
                        {
                            if (!ingroup) _addpermutation(ring,workperm,n);
                            else  _addpermutationunmarked(ring,workperm,n);
                            ingroup = true;
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
        changed = true;
        _addpermutation(ring,p,n);
    }

    return changed;
}

/************************************************************************/

bool _expandschreier(sequential_schreier *gp, sequential_permnode **ring, int n, int* minimal_base) {
    int i,j,nfails,wordlen,skips;
    bool changed;
    sequential_permnode *pn;
#if !MAXN
    DYNALLOC1(int,workperm2,workperm2_sz,n,"expandschreier");
#endif

    pn = *ring;
    if (pn == NULL) return false;

    nfails = 0;
    changed = false;

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
        if (_filterschreier(gp,workperm2,ring,true,-1,n,minimal_base))
        {
            changed = true;
            nfails = 0;
        }
        else
            ++nfails;
    }

    return changed;
}

/************************************************************************/

int* _getorbits(int *fix, int nfix, sequential_schreier *gp, sequential_permnode **ring, int n,
                int* minimal_base, int** orbits_sz) {
    int k;
    sequential_schreier *sh,*sha;

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

void _schreier_freedyn(void) {
    DYNFREE(workperm,workperm_sz);
    DYNFREE(workperm2,workperm2_sz);
    DYNFREE(workpermA,workpermA_sz);
    DYNFREE(workpermB,workpermB_sz);
    DYNFREE(workset,workset_sz);
    DYNFREE(workset2,workset2_sz);
    _clearfreelists();
}
