//
// Created by markus on 01.10.19.
//

/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include "my_schreier.h"

long mmultcount = 0;
long mfiltercount = 0;

static permnode id_permnode;
/* represents identity, no actual content, doesn't need TLS_ATTR */
#define ID_PERMNODE (&id_permnode)

DYNALLSTAT(int, workperm, workperm_sz);
DYNALLSTAT(int, workperm2, workperm2_sz);
DYNALLSTAT(int, workpermA, workpermA_sz);
DYNALLSTAT(int, workpermB, workpermB_sz);
DYNALLSTAT(set, workset, workset_sz);
DYNALLSTAT(set, workset2, workset2_sz);

static TLS_ATTR schreier *mschreier_freelist = NULL;
/* Freelist of scheier structures connected by next field.
     * vec, pwr and orbits fields are assumed allocated. */
static TLS_ATTR permnode *mpermnode_freelist = NULL;
/* Freelist of permnode structures connected by next field.
     * p[] is assumed extended. */

static TLS_ATTR int mschreierfails = SCHREIERFAILS;

#define TMP

static boolean mfilterschreier(schreier *, int *, permnode **, boolean, int, int);

#define PNCODE(x) ((int)(((size_t)(x)>>3)&0xFFFUL))

/* #define TESTP(id,p,n) testispermutation(id,p,n) */
#define TESTP(id, p, n)

/************************************************************************/

static void
mtestispermutation(int id, int *p, int n)
/* For debugging purposes, crash with a message if p[0..n-1] is
   not a permutation. */
{
    int i, m;
    DYNALLSTAT(set, seen, seen_sz);

    for (i = 0; i < n; ++i)
        if (p[i] < 0 || p[i] > n) break;

    if (i < n) {
        fprintf(stderr, ">E Bad permutation (id=%d): n=%d p[%d]=%d\n",
                id, n, i, p[i]);
        exit(1);
    }

    m = SETWORDSNEEDED(n);
    DYNALLOC1(set, seen, seen_sz, m, "malloc seen");
    EMPTYSET(seen, m);

    for (i = 0; i < n; ++i) {
        if (ISELEMENT(seen, p[i])) {
            fprintf(stderr,
                    ">E Bad permutation (id=%d): n=%d p[%d]=%d is a repeat\n",
                    id, n, i, p[i]);
            exit(1);
        }
                ADDELEMENT(seen, p[i]);
    }
}

/************************************************************************/

int
mschreier_fails(int nfails)
/* Set the number of consecutive failures for filtering;
 * A value of <= 0 defaults to SCHREIERFAILS.
 * The function value is the previous setting. */
{
    int prev;

    prev = mschreierfails;

    if (nfails <= 0) mschreierfails = SCHREIERFAILS;
    else mschreierfails = nfails;

    return prev;
}

/************************************************************************/

static void
mclearfreelists(void)
/* Clear the schreier and permnode freelists */
{
    schreier *sh, *nextsh;
    permnode *p, *nextp;

    nextsh = mschreier_freelist;
    while (nextsh) {
        sh = nextsh;
        nextsh = sh->next;
        free(sh->vec);
        free(sh->pwr);
        free(sh->orbits);
        free(sh);
    }
    mschreier_freelist = NULL;

    nextp = mpermnode_freelist;
    while (nextp) {
        p = nextp;
        nextp = p->next;
        free(p);
    }
    mpermnode_freelist = NULL;
}

/************************************************************************/

static permnode
*mnewpermnode(int n)
/* Allocate a new permode structure, with initialized next/prev fields */
{
    permnode *p;

    while (mpermnode_freelist) {
        p = mpermnode_freelist;
        mpermnode_freelist = p->next;
        if (p->nalloc >= n && p->nalloc <= n + 100) {
            p->next = p->prev = NULL;
            p->mark = 0;
            return p;
        } else
            free(p);
    }

    p = (permnode *) malloc(sizeof(permnode) + (n - 2) * sizeof(int));

    if (p == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newpermnode()\n");
        exit(1);
    }

    p->next = p->prev = NULL;
    p->nalloc = n;

    return p;
}

/************************************************************************/

static schreier
*mnewschreier(int n)
/* Allocate a new schreier structure, with initialised next field */
{
    schreier *sh;

    while (mschreier_freelist) {
        sh = mschreier_freelist;
        mschreier_freelist = sh->next;
        if (sh->nalloc >= n && sh->nalloc <= n + 100) {
            sh->next = NULL;
            return sh;
        } else {
            free(sh->vec);
            free(sh->pwr);
            free(sh->orbits);
            free(sh);
        }
    }

    sh = (schreier *) malloc(sizeof(schreier));

    if (sh == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->vec = (permnode **) malloc(sizeof(permnode *) * n);
    sh->pwr = (int *) malloc(sizeof(int) * n);
    sh->orbits = (int *) malloc(sizeof(int) * n);

    if (sh->vec == NULL || sh->pwr == NULL || sh->orbits == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->next = NULL;
    sh->nalloc = n;

    return sh;
}

/************************************************************************/

void
mfreeschreier(schreier **gp, permnode **gens)
/* Free schreier structure and permutation ring.  Assume this is everything. */
/* Use NULL for arguments which don't need freeing. */
{
    schreier *sh, *nextsh;
    permnode *p, *nextp;

    if (gp && *gp) {
        nextsh = *gp;
        while (nextsh) {
            sh = nextsh;
            nextsh = sh->next;
            sh->next = mschreier_freelist;
            mschreier_freelist = sh;
        }
        *gp = NULL;
    }

    if (gens && *gens) {
        p = *gens;
        do {
            nextp = p->next;
            p->next = mpermnode_freelist;
            mpermnode_freelist = p;
            p = nextp;
        } while (p != *gens);
        *gens = NULL;
    }
}

/************************************************************************/

permnode *
mfindpermutation(permnode *pn, int *p, int n)
/* Return a pointer to permutation p in the circular list,
 * or NULL if it isn't present. */
{
    permnode *rn;
    int i;

    if (!pn) return NULL;

    rn = pn;
    do {
        for (i = 0; i < n; ++i)
            if (rn->p[i] != p[i]) break;
        if (i == n) return rn;
        rn = rn->next;
    } while (rn != pn);

    return NULL;
}

/************************************************************************/

void
maddpermutation(permnode **ring, int *p, int n)
/* Add new permutation to circular list, marked.
 * and return pointer to it in *ring. */
{
    permnode *pn, *rn;

    pn = mnewpermnode(n);
    rn = *ring;

    memcpy(pn->p, p, n * sizeof(int));

    if (!rn)
        pn->next = pn->prev = pn;
    else {
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
maddpermutationunmarked(permnode **ring, int *p, int n)
/* Add new permutation to circular list, not marked.
 * and return pointer to it in *ring. */
{
    TESTP(3, p, n);
    maddpermutation(ring, p, n);
    (*ring)->mark = 0;
}

/************************************************************************/

boolean
maddgenerator(schreier **gp, permnode **ring, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    TESTP(2, p, n);
    return mfilterschreier(*gp, p, ring, FALSE, -1, n);
}

/************************************************************************/

boolean
mcondaddgenerator(schreier **gp, permnode **ring, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered, but this version will
 * always notice if this permutation precisely is present.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    TESTP(4, p, n);
    if (mfindpermutation(*ring, p, n))
        return FALSE;
    else
        return mfilterschreier(*gp, p, ring, FALSE, -1, n);
}

/************************************************************************/

static void
mdelpermnode(permnode **ring)
/* Delete permnode at head of circular list, making the next node head. */
{
    permnode *newring;

    if (!*ring) return;

    if ((*ring)->next == *ring)
        newring = NULL;
    else {
        newring = (*ring)->next;
        newring->prev = (*ring)->prev;
        (*ring)->prev->next = newring;
    }

    (*ring)->next = mpermnode_freelist;
    mpermnode_freelist = *ring;

    *ring = newring;
}

/************************************************************************/

void
mdeleteunmarked(permnode **ring)
/* Delete all permutations in the ring that are not marked */
{
    permnode *pn, *firstmarked;

    pn = *ring;
    firstmarked = NULL;

    while (pn != NULL && pn != firstmarked) {
        if (pn->mark) {
            if (!firstmarked) firstmarked = pn;
            pn = pn->next;
        } else
            mdelpermnode(&pn);
    }

    *ring = pn;
}

/************************************************************************/

static void
mclearvector(permnode **vec, permnode **ring, int n)
/* clear vec[0..n-1], freeing permnodes that have no other references
 * and are not marked */
{
    int i;

    for (i = 0; i < n; ++i)
        if (vec[i]) {
            if (vec[i] != ID_PERMNODE) {
                --(vec[i]->refcount);
                if (vec[i]->refcount == 0 && !vec[i]->mark) {
                    *ring = vec[i];
                    mdelpermnode(ring);
                }
            }
            vec[i] = NULL;
        }
}

/************************************************************************/

static void
minitschreier(schreier *sh, int n)
/* Initialise schreier structure to trivial orbits and empty vector */
{
    int i;

    sh->fixed = -1;
    for (i = 0; i < n; ++i) {
        sh->vec[i] = NULL;
        sh->orbits[i] = i;
    }
}

/************************************************************************/

void
mnewgroup(schreier **sh, permnode **ring, int n)
/* Make the trivial group, allow for ring to be set elsewhere */
{
    *sh = mnewschreier(n);
    minitschreier(*sh, n);
    if (ring) *ring = NULL;
}

/************************************************************************/

static void
mapplyperm(int *wp, int *p, int k, int n)
/* Apply the permutation p, k times to each element of wp */
{
    int i, j, cyclen, kk, m;

    TESTP(1, p, n);

    if (k <= 5) {
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
    } else if (k <= 19) {
        DYNALLOC1(int, workpermA, workpermA_sz, n, "applyperm");
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
    } else {
        m = SETWORDSNEEDED(n);
        DYNALLOC1(int, workpermA, workpermA_sz, n, "applyperm");
        DYNALLOC1(int, workpermB, workpermB_sz, n, "applyperm");
        DYNALLOC1(set, workset2, workset2_sz, m, "applyperm");

        EMPTYSET(workset2, m);

        /* We will construct p^k in workpermB one cycle at a time. */

        for (i = 0; i < n; ++i) {
            if (ISELEMENT(workset2, i)) continue;
            if (p[i] == i)
                workpermB[i] = i;
            else {
                cyclen = 1;
                workpermA[0] = i;
                for (j = p[i]; j != i; j = p[j]) {
                    workpermA[cyclen++] = j;
                            ADDELEMENT(workset2, j);
                }
                kk = k % cyclen;
                for (j = 0; j < cyclen; ++j) {
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
mfilterschreier(schreier *gp, int *p, permnode **ring,
                boolean ingroup, int maxlevel, int n)
/* Filter permutation p up to level maxlevel of gp.
 * Use ingroup=TRUE if p is known to be in the group, otherwise
 * at least one equivalent generator is added unless it is proved
 * (nondeterministically) that it is in the group already.
 * maxlevel < 0 means no limit, maxlevel=0 means top level only, etc.
 * Return TRUE iff some change is made. */
{
    int i, j, j1, j2, lev;
    int ipwr;
    schreier *sh;
    int *orbits, *pwr;
    permnode **vec, *curr;
    boolean changed, lchanged, ident;

    DYNALLOC1(int, workperm, workperm_sz, n, "filterschreier");

    ++mfiltercount;

    memcpy(workperm, p, n * sizeof(int));

    if (*ring && p == (*ring)->p) {
        ingroup = TRUE;
        curr = *ring;
    } else
        curr = NULL;

/* curr is the location of workperm in ring, if anywhere */

    sh = gp;
    changed = FALSE;
    if (maxlevel < 0) maxlevel = n + 1;

    for (lev = 0; lev <= maxlevel; ++lev) {
        for (i = 0; i < n; ++i) if (workperm[i] != i) break;
        ident = (i == n);
        if (ident) break;

        lchanged = FALSE;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;
        for (i = 0; i < n; ++i) {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[workperm[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];

            if (j1 != j2) {
                lchanged = TRUE;
                if (j1 < j2) orbits[j2] = j1;
                else orbits[j1] = j2;
            }
        }
        if (lchanged)
            for (i = 0; i < n; ++i) orbits[i] = orbits[orbits[i]];

        if (lchanged) changed = TRUE;

        if (sh->fixed >= 0) {
            for (i = 0; i < n; ++i)
                if (vec[i] && !vec[workperm[i]]) {
                    changed = TRUE;
                    ipwr = 0;
                    for (j = workperm[i]; !vec[j]; j = workperm[j]) ++ipwr;

                    for (j = workperm[i]; !vec[j]; j = workperm[j]) {
                        if (!curr) {
                            if (!ingroup) maddpermutation(ring, workperm, n);
                            else maddpermutationunmarked(ring, workperm, n);
                            ingroup = TRUE;
                            curr = *ring;
                        }
                        vec[j] = curr;
                        pwr[j] = ipwr--;
                        ++curr->refcount;
                    }
                }

            j = workperm[sh->fixed];

            while (j != sh->fixed) {
                mapplyperm(workperm, vec[j]->p, pwr[j], n);
                ++mmultcount;
                curr = NULL;
                j = workperm[sh->fixed];
            }
            sh = sh->next;
        } else
            break;
    }

    if (!ident && !ingroup) {
        changed = TRUE;
        maddpermutation(ring, p, n);
    }

    return changed;
}

/************************************************************************/


static boolean
mfilterschreier_interval(schreier *gp, int *p, permnode **ring,
                boolean ingroup, int maxlevel, int n, int startlevel, int endlevel)
/* Interval version for pipelining:
 * if startlevel = 0,        initialize all DS
 * if endlevel   = maxlevel, delete all DS and return properly
 * if endlevel  != maxlevel, return filter state (workperm, etc.), leave DS initialized
 * */
/* Filter permutation p up to level maxlevel of gp.
 * Use ingroup=TRUE if p is known to be in the group, otherwise
 * at least one equivalent generator is added unless it is proved
 * (nondeterministically) that it is in the group already.
 * maxlevel < 0 means no limit, maxlevel=0 means top level only, etc.
 * Return TRUE iff some change is made. */
{
    int i, j, j1, j2, lev;
    int ipwr;
    schreier *sh;
    int *orbits, *pwr;
    permnode **vec, *curr;
    boolean changed, lchanged, ident;

    DYNALLOC1(int, workperm, workperm_sz, n, "filterschreier");

    ++mfiltercount;

    memcpy(workperm, p, n * sizeof(int));

    if (*ring && p == (*ring)->p) {
        ingroup = TRUE;
        curr = *ring;
    } else
        curr = NULL;

/* curr is the location of workperm in ring, if anywhere */

    sh = gp;
    changed = FALSE;
    if (maxlevel < 0) maxlevel = n + 1;

    for (lev = 0; lev <= maxlevel; ++lev) {
        for (i = 0; i < n; ++i) if (workperm[i] != i) break;
        ident = (i == n);
        if (ident) break;

        lchanged = FALSE;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;
        for (i = 0; i < n; ++i) {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[workperm[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];

            if (j1 != j2) {
                lchanged = TRUE;
                if (j1 < j2) orbits[j2] = j1;
                else orbits[j1] = j2;
            }
        }
        if (lchanged)
            for (i = 0; i < n; ++i) orbits[i] = orbits[orbits[i]];

        if (lchanged) changed = TRUE;

        if (sh->fixed >= 0) {
            for (i = 0; i < n; ++i)
                if (vec[i] && !vec[workperm[i]]) {
                    changed = TRUE;
                    ipwr = 0;
                    for (j = workperm[i]; !vec[j]; j = workperm[j]) ++ipwr;

                    for (j = workperm[i]; !vec[j]; j = workperm[j]) {
                        if (!curr) {
                            if (!ingroup) maddpermutation(ring, workperm, n);
                            else maddpermutationunmarked(ring, workperm, n);
                            ingroup = TRUE;
                            curr = *ring;
                        }
                        vec[j] = curr;
                        pwr[j] = ipwr--;
                        ++curr->refcount;
                    }
                }

            j = workperm[sh->fixed];

            while (j != sh->fixed) {
                mapplyperm(workperm, vec[j]->p, pwr[j], n);
                ++mmultcount;
                curr = NULL;
                j = workperm[sh->fixed];
            }
            sh = sh->next;
        } else
            break;
    }

    if (!ident && !ingroup) {
        changed = TRUE;
        maddpermutation(ring, p, n);
    }

    return changed;
}

/************************************************************************/

boolean
mexpandschreier(schreier *gp, permnode **ring, int n)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i, j, nfails, wordlen, skips;
    boolean changed;
    permnode *pn;

    DYNALLOC1(int, workperm2, workperm2_sz, n, "expandschreier");

    pn = *ring;
    if (pn == NULL) return FALSE;

    nfails = 0;
    changed = FALSE;

    for (skips = KRAN(17); --skips >= 0;) pn = pn->next;

    memcpy(workperm2, pn->p, n * sizeof(int));

    while (nfails < mschreierfails) {
        wordlen = 1 + KRAN(3);
        for (j = 0; j < wordlen; ++j) {
            for (skips = KRAN(17); --skips >= 0;) pn = pn->next;
            for (i = 0; i < n; ++i) workperm2[i] = pn->p[workperm2[i]];
        }
        if (mfilterschreier(gp, workperm2, ring, TRUE, -1, n)) {
            changed = TRUE;
            nfails = 0;
        } else
            ++nfails;
    }

    return changed;
}

/************************************************************************/

int *
mgetorbits(int *fix, int nfix, schreier *gp, permnode **ring, int n)
/* Get a pointer to the orbits for this partial base. The pointer
 * remains valid until pruneset(), getorbits(), getorbitsmin()
 * or grouporder() is called with an incompatible base (neither a
 * prefix nor an extension). The contents of the array pointed to
 * MUST NOT BE MODIFIED by the calling program.
 */
{
    int k;
    schreier *sh, *sha;

    sh = gp;
    for (k = 0; k < nfix; ++k) {
        if (sh->fixed != fix[k]) break;
        sh = sh->next;
    }

    if (k == nfix) return sh->orbits;

    sh->fixed = fix[k];
    mclearvector(sh->vec, ring, n);
    sh->vec[fix[k]] = ID_PERMNODE;

    for (sha = sh->next; sha; sha = sha->next) mclearvector(sha->vec, ring, n);

    for (++k; k <= nfix; ++k) {
        if (!sh->next) sh->next = mnewschreier(n);
        sh = sh->next;
        minitschreier(sh, n);
        if (k < nfix) {
            sh->fixed = fix[k];
            sh->vec[fix[k]] = ID_PERMNODE;
        } else
            sh->fixed = -1;
    }

    if (*ring) mexpandschreier(gp, ring, n);
    return sh->orbits;
}

/************************************************************************/

int
mgetorbitsmin(int *fix, int nfix, schreier *gp, permnode **ring,
              int **orbits, int *cell, int ncell, int n, boolean changed)
/* If the basis elements fix[0..nfix-1] are minimal in their orbits,
 * as far as we know, return value nfix and set *orbits to point
 * to orbits fixing fix[0..nfix-1]. If fix[i] is seen to be not
 * minimal for some i <= nfix-1, return i and set *orbits to point
 * to orbits fixing fix[0..i-1]. If the partial base is already
 * known, or fix[0..nfix-1] can already be seen to be non-minimal,
 * do this work without more filtering. This shortcut is turned
 * off if changed==TRUE. Otherwise, filter until schreierfails
 * failures.
 * The pointer returned remains valid until pruneset(), getorbits(),
 * getorbitsmin() or grouporder() is called with an incompatible base
 * (neither a prefix nor an extension). The contents of the array
 * pointed to MUST NOT BE MODIFIED by the calling program.
 * If cell != NULL, return early if possible when cell[0..ncell-1]
 * are all in the same orbit fixing fix[0..nfix-1]. Otherwise
 * cell,ncell play no part in the computation.
 */
{
    schreier *sh, *sha;
    int *fixorbs;
    int i, j, k, icell, nfails, wordlen, skips;
    permnode *pn;

    DYNALLOC1(int, workperm2, workperm2_sz, n, "expandschreier");

    sh = gp;
    k = 0;
    if (!changed)
        for (k = 0; k < nfix; ++k) {
            if (sh->orbits[fix[k]] != fix[k]) {
                *orbits = sh->orbits;
                return k;
            }
            if (sh->fixed != fix[k]) break;
            sh = sh->next;
        }

    if (k == nfix) {
        *orbits = sh->orbits;
        return nfix;
    }

    sh->fixed = fix[k];
    mclearvector(sh->vec, ring, n);
    sh->vec[fix[k]] = ID_PERMNODE;

    for (sha = sh->next; sha; sha = sha->next) mclearvector(sha->vec, ring, n);

    for (++k; k <= nfix; ++k) {
        if (!sh->next) sh->next = mnewschreier(n);
        sh = sh->next;
        minitschreier(sh, n);
        if (k < nfix) {
            sh->fixed = fix[k];
            sh->vec[fix[k]] = ID_PERMNODE;
        } else
            sh->fixed = -1;
    }
    *orbits = fixorbs = sh->orbits;

    if (cell) {
        for (icell = 1; icell < ncell; ++icell)
            if (fixorbs[cell[icell]] != fixorbs[cell[0]]) break;

        if (icell >= ncell) return nfix;
    }

    if (*ring) {
        pn = *ring;

        nfails = 0;

        for (skips = KRAN(17); --skips >= 0;) pn = pn->next;

        memcpy(workperm2, pn->p, n * sizeof(int));

        while (nfails < mschreierfails) {
            wordlen = 1 + KRAN(3);
            for (j = 0; j < wordlen; ++j) {
                for (skips = KRAN(17); --skips >= 0;) pn = pn->next;
                for (i = 0; i < n; ++i) workperm2[i] = pn->p[workperm2[i]];
            }
            if (mfilterschreier(gp, workperm2, ring, TRUE, -1, n)) {
                nfails = 0;
                sh = gp;
                for (k = 0; k < nfix; ++k) {
                    if (sh->orbits[fix[k]] != fix[k]) {
                        *orbits = sh->orbits;
                        return k;
                    }
                    sh = sh->next;
                }
                if (cell) {
                    for (; icell < ncell; ++icell)
                        if (fixorbs[cell[icell]] != fixorbs[cell[0]]) break;

                    if (icell >= ncell) return nfix;
                }
            } else
                ++nfails;
        }
    }

    return nfix;
}

/************************************************************************/

void
mpruneset(set *fixset, schreier *gp, permnode **ring, set *x, int m, int n)
/* Remove from x any point not minimal for the orbits for this base.
 * If the base is already known, just provide the orbits without
 * more filtering. Otherwise, filter until schreierfails failures.
 */
{
    int i, k;
    schreier *sh, *sha;
    int *orbits;

    DYNALLOC1(set, workset, workset_sz, m, "pruneset");
    for (i = 0; i < m; ++i) workset[i] = fixset[i];

    sh = gp;
    while (sh->fixed >= 0 && ISELEMENT(workset, sh->fixed)) {
                DELELEMENT(workset, sh->fixed);
        sh = sh->next;
    }

    k = nextelement(workset, m, -1);
    if (k < 0)
        orbits = sh->orbits;
    else {
        sh->fixed = k;
        mclearvector(sh->vec, ring, n);
        sh->vec[k] = ID_PERMNODE;

        for (sha = sh->next; sha; sha = sha->next)
            mclearvector(sha->vec, ring, n);

        while ((k = nextelement(workset, m, k)) >= 0) {
            if (!sh->next) sh->next = mnewschreier(n);
            sh = sh->next;
            minitschreier(sh, n);
            sh->fixed = k;
            sh->vec[k] = ID_PERMNODE;
        }
        if (!sh->next) sh->next = mnewschreier(n);
        sh = sh->next;
        minitschreier(sh, n);
        sh->fixed = -1;

        if (*ring) mexpandschreier(gp, ring, n);
        orbits = sh->orbits;
    }

    for (k = -1; (k = nextelement(x, m, k)) >= 0;)
        if (orbits[k] != k) DELELEMENT (x, k);
}

/************************************************************************/
int
mschreier_gens(permnode *ring)
/* Returns the number of generators in the ring */
{
    int j;
    permnode *pn;

    if (!ring) j = 0;
    else for (j = 1, pn = ring->next; pn != ring; pn = pn->next) ++j;

    return j;
}

void
mgrouporder(int *fix, int nfix, schreier *gp, permnode **ring,
            double *grpsize1, int *grpsize2, int n)
/* process the base like in getorbits(), then return the product of the
 * orbits along the base, using the largest orbit at the end if the
 * base is not complete.
*/
{
    schreier *sh;
    int i, j, k, fx;
    int *orb;

    DYNALLOC1(int, workperm, workperm_sz, n, "grouporder");

    mgetorbits(fix, nfix, gp, ring, n);
    mexpandschreier(gp, ring, n);
    mexpandschreier(gp, ring, n);
    *grpsize1 = 1.0;
    *grpsize2 = 0;

    for (i = 0, sh = gp; i < nfix; ++i, sh = sh->next) {
        orb = sh->orbits;
        fx = orb[sh->fixed];
        k = 0;
        for (j = fx; j < n; ++j) if (orb[j] == fx) ++k;
        MULTIPLY(*grpsize1, *grpsize2, k);
    }

    orb = sh->orbits;
    k = 1;
    for (i = 0; i < n; ++i)
        if (orb[i] == i)
            workperm[i] = 1;
        else {
            ++workperm[orb[i]];
            if (workperm[orb[i]] > k) k = workperm[orb[i]];
        }
    MULTIPLY(*grpsize1, *grpsize2, k);
}

/************************************************************************/

void
mschreier_freedyn(void) {
    DYNFREE(workperm, workperm_sz);
    DYNFREE(workperm2, workperm2_sz);
    DYNFREE(workpermA, workpermA_sz);
    DYNFREE(workpermB, workpermB_sz);
    DYNFREE(workset, workset_sz);
    DYNFREE(workset2, workset2_sz);
    mclearfreelists();
}