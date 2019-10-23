//
// Created by markus on 01.10.19.
//

#ifndef DEJAVU_PIPELINE_SCHREIER_H
#define DEJAVU_PIPELINE_SCHREIER_H

#include "nauty/nauty.h"
#include "nauty/naurng.h"

#define DYNALLSTAT_NOSTATIC(type,name,name_sz) \
	type *name; size_t name_sz=0;

typedef struct mpermnodestruct
{
    struct mpermnodestruct *prev,*next;   /* prev&next in circular list */
    unsigned long refcount;              /* number of references */
    int nalloc;                          /* size of p[] in ints,
                                            <= 0 for a perm marker */
    int mark;                            /* a mark, 0 unless changed */
    int p[2];                            /* actual vector, extended to
                                            nalloc enties */
} mpermnode;

typedef struct mschreierlevel
{
    struct mschreierlevel *next;    /* down one level, if any */
    int fixed;                     /* fixed at next level, -1 if none */
    /* Invariant: next=NULL => fixed = -1 */
    int nalloc;                    /* size of vec[] and orbits[] */
    mpermnode **vec;                /* vec[i]^pwr[i] is edge label, */
    int *pwr;                      /*  transitive closure maps i->fixed */
    int *orbits;                   /* vector of orbits */
    int *fixed_orbit = nullptr;
    int  fixed_orbit_sz = -1;
    mpermnode *marker;              /* points to marker for this level */
} mschreier;


typedef struct filterstatestruct {
    int level;
    bool loop_break;
    mschreier *sh;
    int *orbits, *pwr;
    mpermnode **vec, *curr;
    boolean changed, lchanged, ident;
    bool ingroup;
    bool counts_towards_abort;
    int* workperm;
    size_t workperm_sz;
} filterstate;

typedef struct randomelementstruct {
    int* perm;
    size_t perm_sz;
} random_element;

#define SCHREIERFAILS 10
/* Default number of Schreier failures before giving up. */

#ifdef __cplusplus
extern "C" {
#endif

/* See separate file schreier.txt for a description of usage. */

extern void mfreeschreier(mschreier **gp, mpermnode **gens);
extern void maddpermutation(mpermnode **ring, int *p, int n);
extern mpermnode *mfindpermutation(mpermnode *gens, int *p, int n);
extern bool generate_random_element(mschreier *gp, mpermnode **ring, int n, random_element* element);
extern boolean mfilterschreier(mschreier *, int *, mpermnode **, boolean, int, int);
extern boolean mfilterschreier_interval(mschreier *, int *, mpermnode **, boolean, int, int, int, int, filterstate* state);
void free_random_element(random_element* r);
extern boolean maddgenerator(mschreier **gp, mpermnode **gens, int *p, int n);
extern boolean
mcondaddgenerator(mschreier **gp, mpermnode **gens, int *p, int n);
extern boolean mexpandschreier(mschreier *gp, mpermnode **gens, int n);
extern int *mgetorbits(int *fix, int nfix,
                      mschreier *gp, mpermnode **gens, int n);
extern int mgetorbitsmin(int *fix, int nfix, mschreier *gp, mpermnode **gens,
                        int **orbits, int *cell, int ncell, int n, boolean changed);
extern void mpruneset(set *fixset, mschreier *gp, mpermnode **gens,
                     set *x, int m, int n);
extern void mnewgroup(mschreier **gp, mpermnode **gens, int n);
extern void mschreier_freedyn(void);
extern int mschreier_fails(int nfails);
extern void mdumpschreier(FILE *f, mschreier *gp, mpermnode *gens, int n);
extern int mschreier_gens(mpermnode *gens);
extern void mdeleteunmarked(mpermnode **gens);
extern void mgrouporder(int *fix, int nfix,  mschreier *gp, mpermnode **gens,
                       double *grpsize1, int *grpsize2, int n);
extern void mschreier_check(int wordsize, int m, int n, int version);

#ifdef __cplusplus
}
#endif

#endif //DEJAVU_PIPELINE_SCHREIER_H
