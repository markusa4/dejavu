// adapted from schreier.h of nauty 2.5R9

/* schreier.h - Version 1.2 (January 2013) */

#ifndef  _SCHREIER_H_    /* only process this file once */
#define  _SCHREIER_H_

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

typedef struct _schreierlevel
{
    struct _schreierlevel *next;    /* down one level, if any */
    int fixed;                     /* fixed at next level, -1 if none */
                	/* Invariant: next=NULL => fixed = -1 */
    int nalloc;                    /* size of vec[] and orbits[] */
    permnode **vec;                /* vec[i]^pwr[i] is edge label, */
    int *pwr;                      /*  transitive closure maps i->fixed */
    int *orbits;                   /* vector of orbits */
    int *orbits_sz;
    permnode *marker;              /* points to marker for this level */
} _schreier;

#define SCHREIERFAILS 10 
  /* Default number of Schreier failures before giving up. */

#ifdef __cplusplus
extern "C" {
#endif

/* See separate file schreier.txt for a description of usage. */

extern void _freeschreier(_schreier **gp, permnode **gens);
extern void _addpermutation(permnode **ring, int *p, int n);
extern permnode *_findpermutation(permnode *gens, int *p, int n);
extern boolean _addgenerator(_schreier **gp, permnode **gens, int *p, int n);
extern boolean
	_condaddgenerator(_schreier **gp, permnode **gens, int *p, int n);
extern boolean _expandschreier(_schreier *gp, permnode **gens, int n, int* minimal_base);
extern int *_getorbits(int *fix, int nfix,
		 _schreier *gp, permnode **gens, int n, int* minimal_base, int** orbits_sz);
extern void _newgroup(_schreier **sh, permnode **ring, int n);
extern void _schreier_freedyn(void);
extern int _schreier_fails(int nfails);
extern int _schreier_gens(permnode *gens);
extern void _deleteunmarked(permnode **gens);
extern void _grouporder(int *fix, int nfix, _schreier *gp, permnode **gens,
		double *grpsize1, int *grpsize2, int n);
extern void _schreier_check(int wordsize, int m, int n, int version);

#ifdef __cplusplus
}
#endif

#endif  /* _SCHREIER_H_ */
