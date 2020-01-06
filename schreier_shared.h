// adapted from schreier.h of nauty 2.5R9

/* schreier.h - Version 1.2 (January 2013) */

#ifndef DEJAVU_SCHREIER_SHARED_H
#define DEJAVU_SCHREIER_SHARED_H

#include <mutex>
#include <atomic>

#define DYNALL(type,name,name_sz) \
	type *name = NULL; size_t name_sz=0;

enum sift_type {SIFT_UNIFORM, SIFT_NON_UNIFORM, SIFT_RANDOM};

typedef struct shared_permnodestruct {
    struct shared_permnodestruct *prev,*next;
    std::atomic<int> refcount;
    int nalloc;
    int copied;
    int mark;
    int p[2];
} shared_permnode;

typedef struct shared_schreierlevel {
    struct shared_schreierlevel *next;
    int fixed;
    int nalloc;
    shared_permnode **vec;
    int *pwr;
    int *orbits;
    int *fixed_orbit = nullptr;
    int  fixed_orbit_sz = -1;
    std::mutex* level_lock;
    shared_permnode *marker;
} shared_schreier;


typedef struct filterstatestruct {
    int level;
    bool loop_break;
    shared_schreier *sh;
    int *orbits, *pwr;
    shared_permnode **vec, *curr;
    bool changed, lchanged, ident;
    bool ingroup;
    bool counts_towards_abort;
    int* workperm;
    size_t workperm_sz;
    sift_type stype;
} filterstate;

typedef struct randomelementstruct {
    int* perm;
    size_t perm_sz;
} random_element;

#define SCHREIERFAILS 10

void shared_freeschreier(shared_schreier **gp, shared_permnode **gens);
void shared_addpermutation(shared_permnode **ring, int *p, int n);
shared_permnode *shared_findpermutation(shared_permnode *gens, int *p, int n);
bool    generate_random_element(shared_schreier *gp, shared_permnode **ring, int n, random_element* element);
bool mfilterschreier(shared_schreier *, int *, shared_permnode **, bool, int, int);
bool mfilterschreier_interval(shared_schreier *, int *, shared_permnode **, bool, int, int, int, int, filterstate* state);
bool mfilterschreier_shared(shared_schreier *gp, int *p, shared_permnode **ring, bool ingroup, int maxlevel, int n, int startlevel, int endlevel, filterstate* state, int reported_change_level);
void free_random_element(random_element* r);
bool shared_addgenerator(shared_schreier **gp, shared_permnode **gens, int *p, int n);
bool shared_condaddgenerator(shared_schreier **gp, shared_permnode **gens, int *p, int n);
bool shared_expandschreier(shared_schreier *gp, shared_permnode **gens, int n);
int*    shared_getorbits(int *fix, int nfix,
                             shared_schreier *gp, shared_permnode **gens, int n);
void shared_newgroup(shared_schreier **gp, shared_permnode **gens, int n);
void shared_schreier_freedyn(void);
int  shared_schreier_fails(int nfails);
int  shared_schreier_gens(shared_permnode *gens);
void mdeleteunmarked(shared_permnode **gens);

void shared_grouporder(int *fix, int nfix, shared_schreier *gp, shared_permnode **gens,
                              double *grpsize1, int *grpsize2, int n);



#endif //DEJAVU_SCHREIER_SHARED_H
