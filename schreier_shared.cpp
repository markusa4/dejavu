// adapted from schreier.c of nauty 2.5R9

/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include <assert.h>
#include <iostream>
#include <mutex>
#include "schreier_shared.h"

char pad1[128];
std::mutex circ_mutex;
char pad2[128];

static shared_permnode id_permnode;
#define ID_PERMNODE (&id_permnode)
static shared_schreier *mschreier_freelist = NULL;
static shared_permnode *mpermnode_freelist = NULL;
static int mschreierfails = SCHREIERFAILS;

/************************************************************************/

int
shared_schreier_fails(int nfails)
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
/* Clear the schreier and sequential_permnode freelists */
{
    shared_schreier *sh, *nextsh;
    shared_permnode *p, *nextp;

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

static shared_permnode
*mnewpermnode(int n)
/* Allocate a new permode structure, with initialized next fields */
{
    shared_permnode *p;
    p = (shared_permnode *) malloc(sizeof(shared_permnode) + (n - 2) * sizeof(int));

    if (p == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newpermnode()\n");
        exit(1);
    }

    p->next = p->prev = NULL;
    p->nalloc = n;

    return p;
}

/************************************************************************/

static shared_schreier
*mnewschreier(int n)
/* Allocate a new schreier structure, with initialised next field */
{
    shared_schreier *sh;

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

    sh = (shared_schreier *) malloc(sizeof(shared_schreier));

    if (sh == NULL) {
        fprintf(ERRFILE, ">E malloc failed in newschreier()\n");
        exit(1);
    }

    sh->vec = (shared_permnode **) malloc(sizeof(shared_permnode *) * n);
    sh->pwr = (int *) malloc(sizeof(int) * n);
    sh->orbits = (int *) malloc(sizeof(int) * n);

    sh->fixed_orbit = new int[n];
    sh->fixed_orbit_sz = 0;

    sh->level_lock = new std::mutex;

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
shared_freeschreier(shared_schreier **gp, shared_permnode **gens)
/* Free schreier structure and permutation ring.  Assume this is everything. */
/* Use NULL for arguments which don't need freeing. */
{
    shared_schreier *sh, *nextsh;
    shared_permnode *p, *nextp;

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

shared_permnode *
shared_findpermutation(shared_permnode *gens, int *p, int n)
/* Return a pointer to permutation p in the circular list,
 * or NULL if it isn't present. */
{
    shared_permnode *rn;
    int i;

    if (!gens) return NULL;

    rn = gens;
    do {
        for (i = 0; i < n; ++i)
            if (rn->p[i] != p[i]) break;
        if (i == n) return rn;
        rn = rn->next;
    } while (rn != gens);

    return NULL;
}

/************************************************************************/

void
shared_addpermutation(shared_permnode **ring, int *p, int n)
/* Add new permutation to circular list, marked.
 * and return pointer to it in *ring. */
{
    shared_permnode *pn, *rn;

    pn = mnewpermnode(n);
    pn->next = NULL;

    memcpy(pn->p, p, n * sizeof(int));

    rn = *ring;
    if (!rn) {
        pn->next = pn;
    } else {
        pn->next = rn->next;
        rn->next = pn;
    }

    pn->refcount = 0;
    pn->mark = 1;
    pn->copied = 0;
    *ring = pn;
}

/************************************************************************/

static void
maddpermutationunmarked(shared_permnode **ring, int *p, int n)
/* Add new permutation to circular list, not marked.
 * and return pointer to it in *ring. */
{
    shared_addpermutation(ring, p, n);
    (*ring)->mark = 0;
}

/************************************************************************/

bool
shared_addgenerator(shared_schreier **gp, shared_permnode **gens, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    filterstate state;
    mfilterschreier_interval(*gp, p, gens, false, n + 1, n, 0, 10, &state);
    assert(test);
    return mfilterschreier_interval(*gp, p, gens, false, n + 1, n, 11, n + 1, &state);
}

/************************************************************************/

bool
shared_condaddgenerator(shared_schreier **gp, shared_permnode **gens, int *p, int n)
/* Add new permutation to group, unless it is discovered to be
 * already in the group.  It is is possible to be in the group
 * and yet this fact is not discovered, but this version will
 * always notice if this permutation precisely is present.
 * Return TRUE if the generator (or an equivalent) is added or the
 * group knowledge with the current partial base is improved. */
{
    if (shared_findpermutation(*gens, p, n))
        return false;
    else
        return mfilterschreier(*gp, p, gens, false, -1, n);
}

/************************************************************************/

static void
mdelpermnode(shared_permnode **ring)
/* Delete sequential_permnode at head of circular list, making the next node head. */
{
    shared_permnode *newring;

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
mdeleteunmarked(shared_permnode **ring)
/* Delete all permutations in the ring that are not marked */
{
    shared_permnode *pn, *firstmarked;

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
mclearvector(shared_permnode **vec, shared_permnode **ring, int n)
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
minitschreier(shared_schreier *sh, int n)
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
shared_newgroup(shared_schreier **gp, shared_permnode **gens, int n)
/* Make the trivial group, allow for ring to be set elsewhere */
{
    *gp = mnewschreier(n);
    minitschreier(*gp, n);
    if (gens) *gens = NULL;
}

/************************************************************************/

void
mapplyperm(int *wp, int *p, int k, int n)
/* Apply the permutation p, k times to each element of wp */
{
    int i, j, cyclen, kk, m;

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
        DYNALL(int, mworkpermA, mworkpermA_sz);
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
        DYNALL(int, mworkpermA, mworkpermA_sz);
        DYNALL(int, mworkpermB, mworkpermB_sz);
        DYNALL(set, mworkset2, mworkset2_sz);
        DYNALLOC1(int, mworkpermA, mworkpermA_sz, n, "applyperm");
        DYNALLOC1(int, mworkpermB, mworkpermB_sz, n, "applyperm");
        DYNALLOC1(set, mworkset2, mworkset2_sz, m, "applyperm");

        EMPTYSET(mworkset2, m);

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

bool mfilterschreier(shared_schreier *gp, int *p, shared_permnode **ring,
                        bool ingroup, int maxlevel, int n)
/* Filter permutation p up to level maxlevel of gp.
 * Use ingroup=TRUE if p is known to be in the group, otherwise
 * at least one equivalent generator is added unless it is proved
 * (nondeterministically) that it is in the group already.
 * maxlevel < 0 means no limit, maxlevel=0 means top level only, etc.
 * Return TRUE iff some change is made. */
{
    int i, j, j1, j2, lev;
    int ipwr;
    shared_schreier *sh;
    int *orbits, *pwr;
    shared_permnode **vec, *curr;
    bool changed, lchanged, ident;

    DYNALL(int, mworkperm, mworkperm_sz);
    DYNALLOC1(int, mworkperm, mworkperm_sz, n, "filterschreier");

    //++mfiltercount;

    memcpy(mworkperm, p, n * sizeof(int));

    if (*ring && p == (*ring)->p) {
        ingroup = true;
        curr = *ring;
    } else
        curr = NULL;

/* curr is the location of workperm in ring, if anywhere */

    sh = gp;
    changed = false;
    if (maxlevel < 0) maxlevel = n + 1;

    for (lev = 0; lev <= maxlevel; ++lev) {
        for (i = 0; i < n; ++i) if (mworkperm[i] != i) break;
        ident = (i == n);
        if (ident) break;

        lchanged = false;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;
        for (i = 0; i < n; ++i) {
            j1 = orbits[i];
            while (orbits[j1] != j1) j1 = orbits[j1];
            j2 = orbits[mworkperm[i]];
            while (orbits[j2] != j2) j2 = orbits[j2];

            if (j1 != j2) {
                lchanged = true;
                if (j1 < j2) orbits[j2] = j1;
                else orbits[j1] = j2;
            }
        }
        if (lchanged)
            for (i = 0; i < n; ++i) orbits[i] = orbits[orbits[i]];

        if (lchanged) changed = true;

        if (sh->fixed >= 0) {
            for (i = 0; i < n; ++i)
                if (vec[i] && !vec[mworkperm[i]]) {
                    changed = true;
                    ipwr = 0;
                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) ++ipwr;

                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) {
                        if (!curr) {
                            circ_mutex.lock();
                            if (!ingroup) shared_addpermutation(ring, mworkperm, n);
                            else maddpermutationunmarked(ring, mworkperm, n);
                            ingroup = true;
                            curr = *ring;
                            circ_mutex.unlock();
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
        changed = true;
        shared_addpermutation(ring, p, n);
    }

    DYNFREE(mworkperm, mworkperm_sz);
    return changed;
}

/************************************************************************/


bool mfilterschreier_interval(shared_schreier *gp, int *p, shared_permnode **ring,
                                 bool ingroup, int maxlevel, int n, int startlevel, int endlevel, filterstate* state)
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
    shared_schreier *sh;
    int *orbits, *pwr;
    bool loop_break;
    shared_permnode **vec, *curr;
    bool changed, lchanged, ident;

    DYNALL(int, mworkperm, mworkperm_sz);
    if (maxlevel < 0) maxlevel = n + 1;

    if(startlevel == 0) {
        DYNALLOC1(int, mworkperm, mworkperm_sz, n, "filterschreier");
        //++mfiltercount;
        memcpy(mworkperm, p, n * sizeof(int));

        if (*ring && p == (*ring)->p) {
            ingroup = true;
            curr = *ring;
        } else
            curr = NULL;

        sh = gp;
        changed = false;
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

        lchanged = false;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;

        if(lev == 0) { // orbit algorithm
            for (i = 0; i < n; ++i) {
                j1 = orbits[i];
                while (orbits[j1] != j1) j1 = orbits[j1];
                j2 = orbits[mworkperm[i]];
                while (orbits[j2] != j2) j2 = orbits[j2];

                if (j1 != j2) {
                    lchanged = true;
                    if (j1 < j2) orbits[j2] = j1;
                    else orbits[j1] = j2;
                }
            }
            if (lchanged)
                for (i = 0; i < n; ++i) orbits[i] = orbits[orbits[i]];
        }

        if (lchanged) changed = true;

        if (sh->fixed >= 0) {
            if(sh->fixed_orbit_sz == 0) {
                sh->fixed_orbit[0] = sh->fixed;
                sh->fixed_orbit_sz += 1;
            }
            for (int ii = 0; ii < sh->fixed_orbit_sz; ++ii) {// only look at orbit of sh->fixed here instead of entire domain
                int i = sh->fixed_orbit[ii];
                if (vec[i] && !vec[mworkperm[i]]) {
                    // acquire sh->level lock here
                    changed = true;
                    ipwr = 0;
                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) ++ipwr;
                    for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) {
                        circ_mutex.lock();
                        if (!curr) {
                            if (!ingroup) shared_addpermutation(ring, mworkperm, n);
                            else maddpermutationunmarked(ring, mworkperm, n);
                            ingroup = true;
                            curr = *ring;
                        }
                        vec[j] = curr;
                        pwr[j] = ipwr--; // add intermediate perms here since we will reuse them very often? at least on lower levels?
                        ++curr->refcount;
                        circ_mutex.unlock();
                        assert(sh->fixed_orbit_sz < n);
                        assert(sh->fixed_orbit_sz >= 0);
                        sh->fixed_orbit[sh->fixed_orbit_sz] = j;
                        sh->fixed_orbit_sz += 1;
                    }
                }
            }

            j = mworkperm[sh->fixed];
            while (j != sh->fixed) {
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
            changed = true;
            shared_addpermutation(ring, p, n);
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


bool mfilterschreier_shared(shared_schreier *gp, int *p, shared_permnode **ring,
                            bool ingroup, int maxlevel, int n, int startlevel, int endlevel, filterstate* state, int reported_change_level)
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
    shared_schreier *sh;
    int *orbits, *pwr;
    bool loop_break;
    shared_permnode **vec, *curr;
    bool changed, lchanged, ident;
    bool report_changed = false;

    // identity should not allocate memory, hence early out here!
    for (i = 0; i < n; ++i) if (p[i] != i) {break;}
    if(i == n) return false;

    shared_permnode** r = ring;

    DYNALL(int, mworkperm, mworkperm_sz);
    if (maxlevel < 0) maxlevel = n + 1;

    if(startlevel == 0) {
        //mworkperm = p;
        DYNALLOC1(int, mworkperm, mworkperm_sz, n, "filterschreier");
        //++mfiltercount;
        memcpy(mworkperm, p, n * sizeof(int));

        if (*ring && p == (*ring)->p) {
            ingroup = true;
            curr = *ring;
        } else
            curr = NULL;

        sh = gp;
        changed = false;
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
        if(lev > 0)
            for (i = 0; i < n; ++i) if (mworkperm[i] != i) {break;}
        ident = (i == n);
        if (ident) {loop_break = true;break;}

        lchanged = false;
        orbits = sh->orbits;
        vec = sh->vec;
        pwr = sh->pwr;

        if(lev == 0) { // orbit algorithm
            bool test = sh->level_lock->try_lock();
            if(test) {
                for (i = 0; i < n; ++i) {
                    j1 = orbits[i];
                    while (orbits[j1] != j1) j1 = orbits[j1];
                    assert(i < n && i >= 0);
                    j2 = orbits[mworkperm[i]];
                    while (orbits[j2] != j2) j2 = orbits[j2];

                    if (j1 != j2) {
                        lchanged = true;
                        if (j1 < j2) orbits[j2] = j1;
                        else orbits[j1] = j2;
                    }
                }
                if (lchanged)
                    for (i = 0; i < n; ++i) orbits[i] = orbits[orbits[i]];
                sh->level_lock->unlock();
            }
        }

        if (lchanged) changed = true;

        if (sh->fixed >= 0) {
            if(sh->fixed_orbit_sz == 0) {
                sh->level_lock->lock();
                if(sh->fixed_orbit_sz == 0) {
                    sh->fixed_orbit[0] = sh->fixed;
                    sh->fixed_orbit_sz += 1;
                }
                sh->level_lock->unlock();
            }
            for (int ii = 0; ii < sh->fixed_orbit_sz; ++ii) {
                // only looking at orbit of sh->fixed here instead of entire domain
                int i = sh->fixed_orbit[ii];
                if (vec[i] && !vec[mworkperm[i]]) {
                    // acquire sh->level lock here
                    sh->level_lock->lock();
                    if (vec[i] && !vec[mworkperm[i]]) { // need to check again, though (data race possible)
                        changed = true;
                        ipwr = 0;
                        for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) ++ipwr;
                        for (j = mworkperm[i]; !vec[j]; j = mworkperm[j]) {
                            if (!curr) {
                                circ_mutex.lock();
                                if (!ingroup) shared_addpermutation(r, mworkperm, n);
                                else maddpermutationunmarked(r, mworkperm, n);
                                ingroup = true;
                                curr = *r;
                                circ_mutex.unlock();
                            }
                            pwr[j] = ipwr--;
                            ++curr->refcount;
                            vec[j] = curr;
                            assert(sh->fixed_orbit_sz < n);
                            assert(sh->fixed_orbit_sz >= 0);
                            sh->fixed_orbit[sh->fixed_orbit_sz] = j;
                            sh->fixed_orbit_sz += 1;
                        }
                    }
                    sh->level_lock->unlock();
                }
            }

            j = mworkperm[sh->fixed];
            while (j != sh->fixed) {
                mapplyperm(mworkperm, vec[j]->p, pwr[j], n);
                //++mmultcount;
                curr = NULL;
                j = mworkperm[sh->fixed];
            }
            sh = sh->next;
        } else {
            loop_break = true;
            if(lev <= reported_change_level)
                report_changed = changed;
            break;
        }
        if(lev <= reported_change_level)
            report_changed = changed;
    }

    if(endlevel == maxlevel) {
        if (!ident && !ingroup) {
            std::cout << "should not happen" << std::endl;
            assert(false);
            changed = TRUE;
            report_changed = changed;
            shared_addpermutation(ring, p, n);
        }

        DYNFREE(mworkperm, mworkperm_sz);
        return report_changed;
    } else {
        //assert(false);
        DYNFREE(mworkperm, mworkperm_sz);
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
        return report_changed;
    }
}

/************************************************************************/

bool
shared_expandschreier(shared_schreier *gp, shared_permnode **gens, int n)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i, j, nfails, wordlen, skips;
    bool changed;
    shared_permnode *pn;

    pn = *gens;
    if (pn == NULL) return false;

    DYNALL(int, mworkperm2, mworkperm2_sz);
    DYNALLOC1(int, mworkperm2, mworkperm2_sz, n, "expandschreier");

    nfails = 0;
    changed = false;

    for (skips = KRAN(17); --skips >= 0;) pn = pn->next;

    memcpy(mworkperm2, pn->p, n * sizeof(int));

    while (nfails < mschreierfails) {
        wordlen = 1 + KRAN(3);
        for (j = 0; j < wordlen; ++j) {
            for (skips = KRAN(17); --skips >= 0;) pn = pn->next;
            for (i = 0; i < n; ++i) mworkperm2[i] = pn->p[mworkperm2[i]];
        }
        if (mfilterschreier(gp, mworkperm2, gens, TRUE, -1, n)) {
            changed = TRUE;
            nfails = 0;
        } else
            ++nfails;
    }

    DYNFREE(mworkperm2, mworkperm2_sz);
    return changed;
}

/************************************************************************/

bool generate_random_element(shared_schreier *gp, shared_permnode **ring, int n, random_element* element)
/* filter random elements until schreierfails failures.
 * Return true if it ever expanded. */
{
    int i, j, wordlen, skips;
    shared_permnode *pn;

    circ_mutex.lock();
    pn = *ring;
    if (pn == NULL) {
        circ_mutex.unlock();
        return false;
    }

    DYNALL(int, mworkperm2, mworkperm2_sz);
    mworkperm2 = NULL;
    DYNALLOC1(int, mworkperm2, mworkperm2_sz, n, "expandschreier");

    element->perm    = mworkperm2;
    element->perm_sz = mworkperm2_sz;

    for (skips = KRAN(17); --skips >= 0;) pn = pn->next;

    memcpy(mworkperm2, pn->p, n * sizeof(int));

    for (int h = 0; h < 10; ++h) {
        wordlen = 1 + KRAN(3);
        for (j = 0; j < wordlen; ++j) {
            for (skips = KRAN(17); --skips >= 0;) pn = pn->next;
            for (i = 0; i < n; ++i) mworkperm2[i] = pn->p[mworkperm2[i]];
        }
    }

    circ_mutex.unlock();
    return true;
}

void free_random_element(random_element* r) {
    DYNFREE(r->perm, r->perm_sz);
}

/************************************************************************/

int *
shared_getorbits(int *fix, int nfix, shared_schreier *gp, shared_permnode **gens, int n)
/* Get a pointer to the orbits for this partial base. The pointer
 * remains valid until pruneset(), getorbits(), getorbitsmin()
 * or grouporder() is called with an incompatible base (neither a
 * prefix nor an extension). The contents of the array pointed to
 * MUST NOT BE MODIFIED by the calling program.
 */
{
    //circ_mutex.lock();
    int k;
    shared_schreier *sh, *sha;

    sh = gp;
    for (k = 0; k < nfix; ++k) {
        if (sh->fixed != fix[k]) break;
        sh = sh->next;
    }

    if (k == nfix)  {
        //circ_mutex.unlock();
        return sh->orbits;
    }

    sh->fixed = fix[k];
    mclearvector(sh->vec, gens, n);
    sh->vec[fix[k]] = ID_PERMNODE;

    for (sha = sh->next; sha; sha = sha->next) mclearvector(sha->vec, gens, n);

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

    if (*gens) shared_expandschreier(gp, gens, n);
    //circ_mutex.unlock();
    return sh->orbits;
}

/************************************************************************/
int
shared_schreier_gens(shared_permnode *gens)
/* Returns the number of generators in the ring */
{
    int j;
    shared_permnode *pn;

    if (!gens) j = 0;
    else for (j = 1, pn = gens->next; pn != gens; pn = pn->next) ++j;

    return j;
}

void
shared_grouporder(int *fix, int nfix, shared_schreier *gp, shared_permnode **gens,
                  double *grpsize1, int *grpsize2, int n)
/* process the base like in getorbits(), then return the product of the
 * orbits along the base, using the largest orbit at the end if the
 * base is not complete.
*/
{
    shared_schreier *sh;
    int i, k;

    *grpsize1 = 1.0;
    *grpsize2 = 0;

    for (i = 0, sh = gp; i < nfix; ++i, sh = sh->next) {
        k = sh->fixed_orbit_sz;
        if(k == 0)
            k = 1;
        MULTIPLY(*grpsize1, *grpsize2, k);
    }
}

/************************************************************************/

void
shared_schreier_freedyn(void) {
    mclearfreelists();
}