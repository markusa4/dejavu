//
// Created by markus on 01.10.19.
//

#ifndef DEJAVU_MY_SCHREIER_H
#define DEJAVU_MY_SCHREIER_H

#include "nauty/nauty.h"
#include "nauty/naurng.h"

typedef struct permnodestruct
{
    struct permnodestruct *prev,*next;   /* prev&next in circular list */
    unsigned long refcount;              /* number of references */
    int nalloc;                          /* size of p[] in ints,
                                            <= 0 for a perm marker */
    int mark;                            /* a mark, 0 unless changed */
    int p[2];                            /* actual vector, extended to
                                            nalloc enties */
} permnode;

typedef struct schreierlevel
{
    struct schreierlevel *next;    /* down one level, if any */
    int fixed;                     /* fixed at next level, -1 if none */
    /* Invariant: next=NULL => fixed = -1 */
    int nalloc;                    /* size of vec[] and orbits[] */
    permnode **vec;                /* vec[i]^pwr[i] is edge label, */
    int *pwr;                      /*  transitive closure maps i->fixed */
    int *orbits;                   /* vector of orbits */
    permnode *marker;              /* points to marker for this level */
} schreier;


typedef struct filterstatestruct {
    int level;
    bool loop_break;
    schreier *sh;
    int *orbits, *pwr;
    permnode **vec, *curr;
    boolean changed, lchanged, ident;
    bool ingroup;
} filterstate;

#define SCHREIERFAILS 10
/* Default number of Schreier failures before giving up. */

#ifdef __cplusplus
extern "C" {
#endif

/* See separate file schreier.txt for a description of usage. */

extern void mfreeschreier(schreier **gp, permnode **gens);
extern void maddpermutation(permnode **ring, int *p, int n);
extern permnode *mfindpermutation(permnode *gens, int *p, int n);
extern boolean maddgenerator(schreier **gp, permnode **gens, int *p, int n);
extern boolean
mcondaddgenerator(schreier **gp, permnode **gens, int *p, int n);
extern boolean mexpandschreier(schreier *gp, permnode **gens, int n);
extern int *mgetorbits(int *fix, int nfix,
                      schreier *gp, permnode **gens, int n);
extern int mgetorbitsmin(int *fix, int nfix, schreier *gp, permnode **gens,
                        int **orbits, int *cell, int ncell, int n, boolean changed);
extern void mpruneset(set *fixset, schreier *gp, permnode **gens,
                     set *x, int m, int n);
extern void mnewgroup(schreier **gp, permnode **gens, int n);
extern void mschreier_freedyn(void);
extern int mschreier_fails(int nfails);
extern void mdumpschreier(FILE *f, schreier *gp, permnode *gens, int n);
extern int mschreier_gens(permnode *gens);
extern void mdeleteunmarked(permnode **gens);
extern void mgrouporder(int *fix, int nfix,  schreier *gp, permnode **gens,
                       double *grpsize1, int *grpsize2, int n);
extern void mschreier_check(int wordsize, int m, int n, int version);

#ifdef __cplusplus
}
#endif

#endif //DEJAVU_MY_SCHREIER_H
