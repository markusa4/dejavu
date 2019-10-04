//
// Created by markus on 01.10.19.
//

/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include <assert.h>
#include <iostream>
#include <mutex>
#include "pipeline_schreier.h"

//long mmultcount = 0;
//long mfiltercount = 0;

std::mutex circ_mutex;

static mpermnode id_permnode;
/* represents identity, no actual content, doesn't need TLS_ATTR */
#define ID_PERMNODE (&id_permnode)

static TLS_ATTR mschreier *mschreier_freelist = NULL;
/* Freelist of scheier structures connected by next field.
     * vec, pwr and orbits fields are assumed allocated. */
static TLS_ATTR mpermnode *mpermnode_freelist = NULL;
/* Freelist of permnode structures connected by next field.
     * p[] is assumed extended. */

static TLS_ATTR int mschreierfails = SCHREIERFAILS;

#define TMP

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
    DYNALLSTAT_NOSTATIC(set, seen, seen_sz);

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
    mschreier *sh, *nextsh;
    mpermnode *p, *nextp;

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

static mpermnode
*mnewpermnode(int n)
/* Allocate a new permode structure, with initialized next/prev fields */
{
    mpermnode *p;

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

    p = (mpermnode *) malloc(sizeof(mpermnode) + (n - 2) * sizeof(int));

    if (p == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newpermnode()\n");
        exit(1);
    }

    p->next = p->prev = NULL;
    p->nalloc = n;

    return p;
}

/************************************************************************/

static mschreier
*mnewschreier(int n)
/* Allocate a new schreier structure, with initialised next field */
{
    mschreier *sh;

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

    sh = (mschreier *) malloc(sizeof(mschreier));

    if (sh == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->vec = (mpermnode **) malloc(sizeof(mpermnode *) * n);
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
mfreeschreier(mschreier **gp, mpermnode **gens)
/* Free schreier structure and permutation ring.  Assume this is everything. */
/* Use NULL for arguments which don't need freeing. */
{
    mschreier *sh, *nextsh;
    mpermnode *p, *nextp;

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

mpermnode *
mfindpermutation(mpermnode *pn, int *p, int n)
/* Return a pointer to permutation p in the circular list,
 * or NULL if it isn't present. */
{
    mpermnode *rn;
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
maddpermutation(mpermnode **ring, int *p, int n)
/* Add new permutation to circular list, marked.
 * and return pointer to it in *ring. */
{
    mpermnode *pn, *rn;

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
maddpermutationunmarked(mpermnode **ring, int *p, int n)
/* Add new permutation to circular list, not marked.
 * and return pointer to it in *ring. */
{
    TESTP(3, p, n);
    maddpermutation(ring, p, n);
    (*ring)->mark = 0;
}

/************************************************************************/

boolean
maddgenerator(mschreier **gp, mpermnode **ring, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    filterstate state;
    TESTP(2, p, n);
    bool test = mfilterschreier_interval(*gp, p, ring, FALSE, n + 1, n, 0, 10, &state);
    assert(test);
    return mfilterschreier_interval(*gp, p, ring, FALSE, n + 1, n, 11, n + 1, &state);
}

/************************************************************************/

boolean
mcondaddgenerator(mschreier **gp, mpermnode **ring, int *p, int n)
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
mdelpermnode(mpermnode **ring)
/* Delete permnode at head of circular list, making the next node head. */
{
    mpermnode *newring;

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
mdeleteunmarked(mpermnode **ring)
/* Delete all permutations in the ring that are not marked */
{
    mpermnode *pn, *firstmarked;

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
mclearvector(mpermnode **vec, mpermnode **ring, int n)
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
minitschreier(mschreier *sh, int n)
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
mnewgroup(mschreier **sh, mpermnode **ring, int n)
/* Make the trivial group, allow for ring to be set elsewhere */
{
    *sh = mnewschreier(n);
    minitschreier(*sh, n);
    if (ring) *ring = NULL;
}

/************************************************************************/

void
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
        DYNALLSTAT_NOSTATIC(int, mworkpermA, mworkpermA_sz);
        DYNALLOC1(int, mworkpermA, mworkpermA_sz, n, "applyperm");
        for (i = 0; i < n; ++i) mworkpermA[i] = p[p[p[i]]];
        for (; k >= 6; k -= 6)
            for (i = 0; i < n; ++i) wp[i] = mworkpermA[mworkpermA[wp[i]]];
        if (k == 1)
            for (i = 0; i < n; ++i) wp[i] = p[wp[i]];
        else if (k == 2)
            for (i = 0; i < n; ++i) wp[i] = p[p[wp[i]]];
        else if (k == 3)
            for (i = 0; i < n; ++i) wp[i] = mworkpermA[wp[i]];
        else if (k == 4)
            for (i = 0; i < n; ++i) wp[i] = p[mworkpermA[wp[i]]];
        else if (k == 5)
            for (i = 0; i < n; ++i) wp[i] = p[p[mworkpermA[wp[i]]]];
        DYNFREE(mworkpermA, mworkpermA_sz);
    } else {
        m = SETWORDSNEEDED(n);
        DYNALLSTAT_NOSTATIC(int, mworkpermA, mworkpermA_sz);
        DYNALLSTAT_NOSTATIC(int, mworkpermB, mworkpermB_sz);
        DYNALLSTAT_NOSTATIC(set, mworkset2, mworkset2_sz);
        DYNALLOC1(int, mworkpermA, mworkpermA_sz, n, "applyperm");
        DYNALLOC1(int, mworkpermB, mworkpermB_sz, n, "applyperm");
        DYNALLOC1(set, mworkset2, mworkset2_sz, m, "applyperm");

        EMPTYSET(mworkset2, m); // ToDo: this may be expensive, there should be only 1 workset per thread!

        /* We will construct p^k in workpermB one cycle at a time. */

        for (i = 0; i < n; ++i) {
            if (ISELEMENT(mworkset2, i)) continue;
            if (p[i] == i)
                mworkpermB[i] = i;
            else {
                cyclen = 1;
                mworkpermA[0] = i;
                for (j = p[i]; j != i; j = p[j]) {
                    mworkpermA[cyclen++] = j;
                            ADDELEMENT(mworkset2, j);
                }
                kk = k % cyclen;
                for (j = 0; j < cyclen; ++j) {
                    mworkpermB[mworkpermA[j]] = mworkpermA[kk];
                    if (++kk == cyclen) kk = 0;
                }
            }
        }
        for (i = 0; i < n; ++i) wp[i] = mworkpermB[wp[i]];

        DYNFREE(mworkpermA, mworkpermA_sz);
        DYNFREE(mworkpermB, mworkpermB_sz);
        DYNFREE(mworkset2, mworkset2_sz);
    }
}

/************************************************************************/

boolean mfilterschreier(mschreier *gp, int *p, mpermnode **ring,
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
    mschreier *sh;
    int *orbits, *pwr;
    mpermnode **vec, *curr;
    boolean changed, lchanged, ident;

    DYNALLSTAT_NOSTATIC(int, mworkperm, mworkperm_sz);
    DYNALLOC1(int, mworkperm, mworkperm_sz, n, "filterschreier");

    //++mfiltercount;

    memcpy(mworkperm, p, n * sizeof(int));

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
        for (i = 0; i < n; ++i) if (mworkperm[i] != i) break;
        ident = (i == n);
        if (ident) break;

        lchanged = FALSE;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;
        for (i = 0; i < n; ++i) {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[mworkperm[i]];
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
                if (vec[i] && !vec[mworkperm[i]]) {
                    changed = TRUE;
                    ipwr = 0;
                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) ++ipwr;

                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) {
                        if (!curr) {
                            if (!ingroup) maddpermutation(ring, mworkperm, n);
                            else maddpermutationunmarked(ring, mworkperm, n);
                            ingroup = TRUE;
                            curr = *ring;
                        }
                        vec[j] = curr;
                        pwr[j] = ipwr--;
                        ++curr->refcount;
                    }
                }

            j = mworkperm[sh->fixed];

            while (j != sh->fixed) {
                mapplyperm(mworkperm, vec[j]->p, pwr[j], n);
               // ++mmultcount;
                curr = NULL;
                j = mworkperm[sh->fixed];
            }
            sh = sh->next;
        } else
            break;
    }

    if (!ident && !ingroup) {
        changed = TRUE;
        maddpermutation(ring, p, n);
    }

    DYNFREE(mworkperm, mworkperm_sz);
    return changed;
}

/************************************************************************/


boolean mfilterschreier_interval(mschreier *gp, int *p, mpermnode **ring,
                boolean ingroup, int maxlevel, int n, int startlevel, int endlevel, filterstate* state)
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
    mschreier *sh;
    int *orbits, *pwr;
    bool loop_break;
    mpermnode **vec, *curr;
    boolean changed, lchanged, ident;

    DYNALLSTAT_NOSTATIC(int, mworkperm, mworkperm_sz);
    if (maxlevel < 0) maxlevel = n + 1;

    if(startlevel == 0) {
        DYNALLOC1(int, mworkperm, mworkperm_sz, n, "filterschreier");

        //++mfiltercount;

        memcpy(mworkperm, p, n * sizeof(int));

        if (*ring && p == (*ring)->p) {
            ingroup = TRUE;
            curr = *ring;
        } else
            curr = NULL;

        sh = gp;
        changed = FALSE;
        loop_break = false;
    } else {
        sh = state->sh;
        orbits = state->orbits;
        pwr = state->pwr;
        vec = state->vec;
        curr = state->curr;
        changed = state->changed;
        lchanged = state->lchanged;
        ident = state->ident ;
        assert(state->level == startlevel - 1);
        loop_break = state->loop_break;
        mworkperm = state->workperm;
        mworkperm_sz = state->workperm_sz;
    }

    for (lev = startlevel; lev <= endlevel; ++lev) {
        if(loop_break) {break;}
        for (i = 0; i < n; ++i) if (mworkperm[i] != i) {break;}
        ident = (i == n);
        if (ident) {loop_break = true;break;}

        lchanged = FALSE;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;
        for (i = 0; i < n; ++i) {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[mworkperm[i]];
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
                if (vec[i] && !vec[mworkperm[i]]) {
                    changed = TRUE;
                    ipwr = 0;
                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) ++ipwr;

                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) {
                        circ_mutex.lock();
                        if (!curr) {
                            if (!ingroup) maddpermutation(ring, mworkperm, n);
                            else maddpermutationunmarked(ring, mworkperm, n);
                            ingroup = TRUE;
                            curr = *ring;
                        }
                        vec[j] = curr;
                        pwr[j] = ipwr--;
                        ++curr->refcount;
                        circ_mutex.unlock();
                    }
                }

            j = mworkperm[sh->fixed];
                //int test__ = 0;
            while (j != sh->fixed) {
                /*test__ += 1;
                if(test__ == 1000) {
                    std::cout << "probably infinite loop" << std::endl;
                }*/
                mapplyperm(mworkperm, vec[j]->p, pwr[j], n);
                //++mmultcount;
                curr = NULL;
                j = mworkperm[sh->fixed];
            }
            sh = sh->next;
        } else {
            loop_break = true;
            break;
        }
    }

    if(endlevel == maxlevel) {
        if (!ident && !ingroup) {
            changed = TRUE;
            maddpermutation(ring, p, n);
        }

        DYNFREE(mworkperm, mworkperm_sz);
        return changed;
    } else {
        //assert(false);
        state->sh   = sh;
        state->orbits = orbits;
        state->pwr  = pwr;
        state->vec  = vec;
        state->curr = curr;
        state->changed  = changed;
        state->lchanged = lchanged;
        state->ident    = ident;
        state->level = endlevel;
        state->loop_break = loop_break;
        state->workperm = mworkperm;
        state->workperm_sz = mworkperm_sz;
        return true;
    }
}

/************************************************************************/

boolean
mexpandschreier(mschreier *gp, mpermnode **ring, int n)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i, j, nfails, wordlen, skips;
    boolean changed;
    mpermnode *pn;

    pn = *ring;
    if (pn == NULL) return FALSE;

    DYNALLSTAT_NOSTATIC(int, mworkperm2, mworkperm2_sz);
    DYNALLOC1(int, mworkperm2, mworkperm2_sz, n, "expandschreier");

    nfails = 0;
    changed = FALSE;

    for (skips = KRAN(17); --skips >= 0;) pn = pn->next;

    memcpy(mworkperm2, pn->p, n * sizeof(int));

    while (nfails < mschreierfails) {
        wordlen = 1 + KRAN(3);
        for (j = 0; j < wordlen; ++j) {
            for (skips = KRAN(17); --skips >= 0;) pn = pn->next;
            for (i = 0; i < n; ++i) mworkperm2[i] = pn->p[mworkperm2[i]];
        }
        if (mfilterschreier(gp, mworkperm2, ring, TRUE, -1, n)) {
            changed = TRUE;
            nfails = 0;
        } else
            ++nfails;
    }

    DYNFREE(mworkperm2, mworkperm2_sz);
    return changed;
}

/************************************************************************/

bool generate_random_element(mschreier *gp, mpermnode **ring, int n, random_element* element)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i, j, wordlen, skips;
    boolean changed;
    mpermnode *pn;

    circ_mutex.lock();
    pn = *ring;
    if (pn == NULL) {
        circ_mutex.unlock();
        return false;
    }

    DYNALLSTAT_NOSTATIC(int, mworkperm2, mworkperm2_sz);
    DYNALLOC1(int, mworkperm2, mworkperm2_sz, n, "expandschreier");

    element->perm    = mworkperm2;
    element->perm_sz = mworkperm2_sz;

    changed = FALSE;

    for (skips = KRAN(17); --skips >= 0;) pn = pn->next;

    memcpy(mworkperm2, pn->p, n * sizeof(int));


    wordlen = 1 + KRAN(3);
    for (j = 0; j < wordlen; ++j) {
        for (skips = KRAN(17); --skips >= 0;) pn = pn->next;
        for (i = 0; i < n; ++i) mworkperm2[i] = pn->p[mworkperm2[i]];
    }

    circ_mutex.unlock();
    return true;
}

void free_random_element(random_element* r) {
    DYNFREE(r->perm, r->perm_sz);
}

/************************************************************************/

int *
mgetorbits(int *fix, int nfix, mschreier *gp, mpermnode **ring, int n)
/* Get a pointer to the orbits for this partial base. The pointer
 * remains valid until pruneset(), getorbits(), getorbitsmin()
 * or grouporder() is called with an incompatible base (neither a
 * prefix nor an extension). The contents of the array pointed to
 * MUST NOT BE MODIFIED by the calling program.
 */
{
    circ_mutex.lock();
    int k;
    mschreier *sh, *sha;

    sh = gp;
    for (k = 0; k < nfix; ++k) {
        if (sh->fixed != fix[k]) break;
        sh = sh->next;
    }

    if (k == nfix)  {
        circ_mutex.unlock();
        return sh->orbits;
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

    if (*ring) mexpandschreier(gp, ring, n);
    circ_mutex.unlock();
    return sh->orbits;
}

/************************************************************************/
int
mschreier_gens(mpermnode *ring)
/* Returns the number of generators in the ring */
{
    int j;
    mpermnode *pn;

    if (!ring) j = 0;
    else for (j = 1, pn = ring->next; pn != ring; pn = pn->next) ++j;

    return j;
}

void
mgrouporder(int *fix, int nfix, mschreier *gp, mpermnode **ring,
            double *grpsize1, int *grpsize2, int n)
/* process the base like in getorbits(), then return the product of the
 * orbits along the base, using the largest orbit at the end if the
 * base is not complete.
*/
{
    mschreier *sh;
    int i, j, k, fx;
    int *orb;

    DYNALLSTAT_NOSTATIC(int, mworkperm, mworkperm_sz);
    DYNALLOC1(int, mworkperm, mworkperm_sz, n, "grouporder");

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
            mworkperm[i] = 1;
        else {
            ++mworkperm[orb[i]];
            if (mworkperm[orb[i]] > k) k = mworkperm[orb[i]];
        }
    MULTIPLY(*grpsize1, *grpsize2, k);
    DYNFREE(mworkperm, mworkperm_sz);
}

/************************************************************************/

void
mschreier_freedyn(void) {
    mclearfreelists();
}