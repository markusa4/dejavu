// adapted from schreier.h of nauty 2.5R9

/* schreier.h - Version 1.2 (January 2013) */

#ifndef DEJAVU_SCHREIER_SEQUENTIAL_H
#define DEJAVU_SCHREIER_SEQUENTIAL_H

#include "nauty/nauty.h"
#include "nauty/naurng.h"

typedef struct sequential_permnodestruct {
    struct sequential_permnodestruct *prev,*next;
    unsigned long refcount;
    int nalloc;
    int mark;
    int p[2];
} sequential_permnode;

typedef struct sequential_schreierlevel {
    struct sequential_schreierlevel *next;
    int fixed;
    int nalloc;
    sequential_permnode **vec;
    int *pwr;
    int *orbits;
    int *orbits_sz;
    sequential_permnode *marker;
} sequential_schreier;

#define SCHREIERFAILS 10

extern void _freeschreier(sequential_schreier **gp, sequential_permnode **gens);
extern void _addpermutation(sequential_permnode **ring, int *p, int n);
extern sequential_permnode *_findpermutation(sequential_permnode *gens, int *p, int n);
extern bool _addgenerator(sequential_schreier **gp, sequential_permnode **gens, int *p, int n);
extern bool
	_condaddgenerator(sequential_schreier **gp, sequential_permnode **gens, int *p, int n);
extern bool _expandschreier(sequential_schreier *gp, sequential_permnode **gens, int n, int* minimal_base);
extern int *_getorbits(int *fix, int nfix,
                       sequential_schreier *gp, sequential_permnode **gens, int n, int* minimal_base, int** orbits_sz);
extern void _newgroup(sequential_schreier **sh, sequential_permnode **ring, int n);
extern void _schreier_freedyn(void);
extern int _schreier_fails(int nfails);
extern int _schreier_gens(sequential_permnode *gens);
extern void _deleteunmarked(sequential_permnode **gens);
extern void _grouporder(int *fix, int nfix, sequential_schreier *gp, sequential_permnode **gens,
                        double *grpsize1, int *grpsize2, int n);
extern void _schreier_check(int wordsize, int m, int n, int version);

#endif  /* DEJAVU_SCHREIER_SEQUENTIAL_H */
