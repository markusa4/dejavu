// adapted from schreier.c of nauty 2.5R9

/* schreier.c - procedures for manipulating a permutation group using
 * the random schreier algorithm.  There is a separate file schreier.txt
 * which describes the usage.
 *
 * Written for nauty and traces, Brendan McKay 2010-2013.
 */

#include <assert.h>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <mutex>
#include "schreier_shared.h"
#include "naurng.h"
#include "naudefs.h"

char pad1[128];
std::mutex circ_mutex;
char pad2[128];

static shared_permnode id_permnode;

static shared_schreier *mschreier_freelist = nullptr;
static shared_permnode *mpermnode_freelist = nullptr;
static int mschreierfails = SCHREIERFAILS;

/************************************************************************/

int shared_schreier_fails(int nfails) {
    int prev;

    prev = mschreierfails;

    if (nfails <= 0) mschreierfails = SCHREIERFAILS;
    else mschreierfails = nfails;

    return prev;
}

/************************************************************************/

static void mclearfreelists(void) {
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
    mschreier_freelist = nullptr;

    nextp = mpermnode_freelist;
    while (nextp) {
        p = nextp;
        nextp = p->next;
        free(p);
    }
    mpermnode_freelist = nullptr;
}

/************************************************************************/

static shared_permnode *mnewpermnode(int n) {
    shared_permnode *p;
    p = (shared_permnode *) malloc(sizeof(shared_permnode) + (n - 2) * sizeof(int));

    if (p == nullptr) {
        exit(1);
    }

    p->next = p->prev = nullptr;
    p->nalloc = n;

    return p;
}

/************************************************************************/

static shared_schreier *mnewschreier(int n) {
    shared_schreier *sh;

    while (mschreier_freelist) {
        sh = mschreier_freelist;
        mschreier_freelist = sh->next;
        if (sh->nalloc >= n && sh->nalloc <= n + 100) {
            sh->next = nullptr;
            return sh;
        } else {
            free(sh->vec);
            free(sh->pwr);
            free(sh->orbits);
            free(sh);
        }
    }

    sh = (shared_schreier *) malloc(sizeof(shared_schreier));

    if (sh == nullptr) {
        exit(1);
    }

    sh->vec = (shared_permnode **) malloc(sizeof(shared_permnode *) * n);
    sh->pwr = (int *) malloc(sizeof(int) * n);
    sh->orbits = (int *) malloc(sizeof(int) * n);

    sh->fixed_orbit = new int[n];
    sh->fixed_orbit_sz = 0;

    sh->level_lock = new std::mutex;

    if (sh->vec == nullptr || sh->pwr == nullptr || sh->orbits == nullptr) {
        exit(1);
    }

    sh->next = nullptr;
    sh->nalloc = n;

    return sh;
}

/************************************************************************/

void shared_freeschreier(shared_schreier **gp, shared_permnode **gens) {
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
        *gp = nullptr;
    }

    if (gens && *gens) {
        p = *gens;
        do {
            nextp = p->next;
            p->next = mpermnode_freelist;
            mpermnode_freelist = p;
            p = nextp;
        } while (p != *gens);
        *gens = nullptr;
    }
}

/************************************************************************/

shared_permnode* shared_findpermutation(shared_permnode *gens, int *p, int n) {
    shared_permnode *rn;
    int i;

    if (!gens) return nullptr;

    rn = gens;
    do {
        for (i = 0; i < n; ++i)
            if (rn->p[i] != p[i]) break;
        if (i == n) return rn;
        rn = rn->next;
    } while (rn != gens);

    return nullptr;
}

/************************************************************************/

void shared_addpermutation(shared_permnode **ring, int *p, int n) {
    shared_permnode *pn, *rn;

    pn = mnewpermnode(n);
    pn->next = nullptr;

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

static void maddpermutationunmarked(shared_permnode **ring, int *p, int n) {
    shared_addpermutation(ring, p, n);
    (*ring)->mark = 0;
}

/************************************************************************/

bool shared_addgenerator(shared_schreier **gp, shared_permnode **gens, int *p, int n) {
    filterstate state;
    mfilterschreier_interval(*gp, p, gens, false, n + 1, n, 0, 10, &state);
    return mfilterschreier_interval(*gp, p, gens, false, n + 1, n, 11, n + 1, &state);
}

/************************************************************************/

bool shared_condaddgenerator(shared_schreier **gp, shared_permnode **gens, int *p, int n) {
    if (shared_findpermutation(*gens, p, n))
        return false;
    else
        return mfilterschreier(*gp, p, gens, false, -1, n);
}

/************************************************************************/

static void mdelpermnode(shared_permnode **ring) {
    shared_permnode *newring;

    if (!*ring) return;

    if ((*ring)->next == *ring)
        newring = nullptr;
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

void mdeleteunmarked(shared_permnode **ring) {
    shared_permnode *pn, *firstmarked;

    pn = *ring;
    firstmarked = nullptr;

    while (pn != nullptr && pn != firstmarked) {
        if (pn->mark) {
            if (!firstmarked) firstmarked = pn;
            pn = pn->next;
        } else
            mdelpermnode(&pn);
    }

    *ring = pn;
}

/************************************************************************/

static void mclearvector(shared_permnode **vec, shared_permnode **ring, int n) {
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
            vec[i] = nullptr;
        }
}

/************************************************************************/

static void minitschreier(shared_schreier *sh, int n) {
    int i;

    sh->fixed = -1;
    for (i = 0; i < n; ++i) {
        sh->vec[i] = nullptr;
        sh->orbits[i] = i;
    }
}

/************************************************************************/

void shared_newgroup(shared_schreier **gp, shared_permnode **gens, int n) {
    *gp = mnewschreier(n);
    minitschreier(*gp, n);
    if (gens) *gens = nullptr;
}

/************************************************************************/

void mapplyperm(int *wp, int *p, int k, int n) {
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

bool mfilterschreier(shared_schreier *gp, int *p, shared_permnode **ring, bool ingroup, int maxlevel, int n) {
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
        curr = nullptr;

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
                curr = nullptr;
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


bool mfilterschreier_interval(shared_schreier *gp, int *p, shared_permnode **ring, bool ingroup, int maxlevel,
                              int n, int startlevel, int endlevel, filterstate* state) {
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
            curr = nullptr;

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
            for (int ii = 0; ii < sh->fixed_orbit_sz; ++ii) {
                // only look at orbit of sh->fixed here instead of entire domain
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
                        pwr[j] = ipwr--;
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
                curr = nullptr;
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


// Interval version for pipelining:
// if startlevel = 0,        initialize all DS
// if endlevel   = maxlevel, delete all DS and return properly
// if endlevel  != maxlevel, return filter state (workperm, etc.), leave DS initialized
bool mfilterschreier_shared(shared_schreier *gp, int *p, shared_permnode **ring, bool ingroup, int maxlevel,
                            int n, int startlevel, int endlevel, filterstate* state, int reported_change_level) {
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
        mworkperm = p;
        //DYNALLOC1(int, mworkperm, mworkperm_sz, n, "filterschreier");
        //++mfiltercount;
        //memcpy(mworkperm, p, n * sizeof(int));

        if (*ring && p == (*ring)->p) {
            ingroup = true;
            curr = *ring;
        } else
            curr = nullptr;

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
                curr = nullptr;
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
            changed = true;
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

bool shared_expandschreier(shared_schreier *gp, shared_permnode **gens, int n) {
    int i, j, nfails, wordlen, skips;
    bool changed;
    shared_permnode *pn;

    pn = *gens;
    if (pn == nullptr) return false;

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
        if (mfilterschreier(gp, mworkperm2, gens, true, -1, n)) {
            changed = true;
            nfails = 0;
        } else
            ++nfails;
    }

    DYNFREE(mworkperm2, mworkperm2_sz);
    return changed;
}

/************************************************************************/

bool generate_random_element(shared_schreier *gp, shared_permnode **ring, int n, random_element* element) {
    int i, j, wordlen, skips;
    shared_permnode *pn;

    circ_mutex.lock();
    pn = *ring;
    if (pn == nullptr) {
        circ_mutex.unlock();
        return false;
    }

    DYNALL(int, mworkperm2, mworkperm2_sz);
    mworkperm2 = nullptr;
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

int * shared_getorbits(int *fix, int nfix, shared_schreier *gp, shared_permnode **gens, int n) {
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
int shared_schreier_gens(shared_permnode *gens) {
    int j;
    shared_permnode *pn;

    if (!gens) j = 0;
    else for (j = 1, pn = gens->next; pn != gens; pn = pn->next) ++j;

    return j;
}

void shared_grouporder(int *fix, int nfix, shared_schreier *gp, shared_permnode **gens, double *grpsize1,
                       int *grpsize2, int n) {
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