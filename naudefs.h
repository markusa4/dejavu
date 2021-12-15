// extracted from "nauty.h" of nauty 2.5R9

#ifndef DEJAVU_NAUDEFS_H
#define DEJAVU_NAUDEFS_H

#include <sys/types.h>
//#include <unistd.h>
#include <stddef.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define ALLOCS(x,y) malloc((size_t)(x)*(size_t)(y))
#define REALLOCS(p,x) realloc(p,(size_t)(x))
#define FREES(p) free(p)
#define ID_PERMNODE (&id_permnode)
#define DYNFREE(name,name_sz) \
  { if (name) FREES(name); name = NULL; name_sz = 0;}
#define DYNALLOC1(type,name,name_sz,sz,msg) \
 if ((size_t)(sz) > name_sz) \
 { if (name_sz) FREES(name); name_sz = (sz); \
 if ((name=(type*)ALLOCS(sz,sizeof(type))) == NULL) {}}
#define DYNALLSTAT(type,name,name_sz) \
	static type *name; static size_t name_sz=0
#define MULTIPLY(s1,s2,i) if ((s1 *= i) >= 1e10) {s1 /= 1e10; s2 += 10;}

#endif //DEJAVU_NAUDEFS_H
