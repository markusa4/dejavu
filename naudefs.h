// extracted from "nauty.h" of nauty distribution

#ifndef DEJAVU_NAUDEFS_H
#define DEJAVU_NAUDEFS_H

#include <sys/types.h>
#include <unistd.h>
#include <stddef.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define SETWD(pos) ((pos)>>6)
#define SETBT(pos) ((pos)&0x3F)
#define TIMESWORDSIZE(w) ((w)<<6)    /* w*WORDSIZE */
#define SETWORDSNEEDED(n) ((((n)-1)>>6)+1)
#define BITT bit
typedef unsigned long setword;
typedef setword set;

static const setword bit[] = {01000000000000000000000,0400000000000000000000,
                              0200000000000000000000,0100000000000000000000,
                              040000000000000000000,020000000000000000000,
                              010000000000000000000,04000000000000000000,
                              02000000000000000000,01000000000000000000,
                              0400000000000000000,0200000000000000000,
                              0100000000000000000,040000000000000000,020000000000000000,
                              010000000000000000,04000000000000000,02000000000000000,
                              01000000000000000,0400000000000000,0200000000000000,
                              0100000000000000,040000000000000,020000000000000,
                              010000000000000,04000000000000,02000000000000,
                              01000000000000,0400000000000,0200000000000,0100000000000,
                              040000000000,020000000000,010000000000,04000000000,
                              02000000000,01000000000,0400000000,0200000000,0100000000,
                              040000000,020000000,010000000,04000000,02000000,01000000,
                              0400000,0200000,0100000,040000,020000,010000,04000,
                              02000,01000,0400,0200,0100,040,020,010,04,02,01};

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
#define ADDELEMENT0(setadd,pos)  ((setadd)[SETWD(pos)] |= BITT[SETBT(pos)])
#define DELELEMENT0(setadd,pos)  ((setadd)[SETWD(pos)] &= ~BITT[SETBT(pos)])
#define FLIPELEMENT0(setadd,pos) ((setadd)[SETWD(pos)] ^= BITT[SETBT(pos)])
#define ISELEMENT0(setadd,pos) (((setadd)[SETWD(pos)] & BITT[SETBT(pos)]) != 0)
#define EMPTYSET0(setadd,m) \
    {setword *es; \
    for (es = (setword*)(setadd)+(m); --es >= (setword*)(setadd);) *es=0;}
#define ADDELEMENT ADDELEMENT0
#define DELELEMENT DELELEMENT0
#define FLIPELEMENT FLIPELEMENT0
#define ISELEMENT ISELEMENT0
#define EMPTYSET EMPTYSET0

#define MULTIPLY(s1,s2,i) if ((s1 *= i) >= 1e10) {s1 /= 1e10; s2 += 10;}

#endif //DEJAVU_NAUDEFS_H
