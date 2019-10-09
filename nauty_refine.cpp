//
// Created by markus on 09.10.19.
//

#include "nauty_refine.h"

#define TMP

/*   #define ONE_WORD_SETS  not sure about this!  See notes.txt.  */

/* macros for hash-codes: */
#define MASH(l,i) ((((l) ^ 065435) + (i)) & 077777)
/* : expression whose long value depends only on long l and int/long i.
     Anything goes, preferably non-commutative. */

#define CLEANUP(l) ((int)((l) % 077777))
/* : expression whose value depends on long l and is less than 077777
     when converted to int then short.  Anything goes. */

#if  MAXM==1
#define M 1
#else
#define M m
#endif

#define ACCUM(x,y)   x = (((x) + (y)) & 077777)

static const int fuzz1[] = {037541,061532,005257,026416};
static const int fuzz2[] = {006532,070236,035523,062437};

#define FUZZ1(x) ((x) ^ fuzz1[(x)&3])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])

/* aproto: header new_nauty_protos.h */

//dispatchvec dispatch_sparse =
//        {isautom_sg,testcanlab_sg,updatecan_sg,refine_sg,refine_sg,cheapautom_sg,
//         targetcell_sg,nausparse_freedyn,nausparse_check,init_sg,NULL};

#if !MAXN
#else
static TLS_ATTR short vmark1[MAXN];
static TLS_ATTR short vmark2[MAXN];
static TLS_ATTR int work1[MAXN];
static TLS_ATTR int work2[MAXN];
static TLS_ATTR int work3[MAXN];
static TLS_ATTR int work4[MAXN];
static TLS_ATTR set snwork[40*MAXM];
#endif

#define MARK1(i) W->vmark1[i] = W->vmark1_val
#define UNMARK1(i) vmark1[i] = 0
#define ISMARKED1(i) (vmark1[i] == vmark1_val)
#define ISNOTMARKED1(i) (W->vmark1[i] != W->vmark1_val)

#define MARK2(i) W->vmark2[i] = W->vmark2_val
#define UNMARK2(i) vmark2[i] = 0
#define ISMARKED2(i) (W->vmark2[i] == W->vmark2_val)
#define ISNOTMARKED2(i) (vmark2[i] != vmark2_val)

#if !MAXN
#define RESETMARKS1 {if (W->vmark1_val++ >= 32000) \
    {size_t ij; for (ij=0;ij<W->vmark1_sz;++ij) W->vmark1[ij]=0; W->vmark1_val=1;}}
#define PREPAREMARKS1(nn, vmark1, vmark1_sz, vmark1_val) preparemarks1(nn, vmark1, vmark1_sz, vmark1_val)
#define RESETMARKS2 {if (W->vmark2_val++ >= 32000) \
    {size_t ij; for (ij=0;ij<W->vmark2_sz;++ij) W->vmark2[ij]=0; W->vmark2_val=1;}}
#define PREPAREMARKS2(nn, vmark2, vmark2_sz, vmark2_val) preparemarks2(nn, vmark2, vmark2_sz, vmark2_val)
#else
#define RESETMARKS1 {if (vmark1_val++ >= 32000) \
    {size_t ij; for (ij=0;ij<MAXN;++ij) vmark1[ij]=0; vmark1_val=1;}}
#define PREPAREMARKS1(nn)
#define RESETMARKS2 {if (vmark2_val++ >= 32000) \
    {size_t ij; for (ij=0;ij<MAXN;++ij) vmark2[ij]=0; vmark2_val=1;}}
#define PREPAREMARKS2(nn)
#endif


#define SORT_OF_SORT 3
#define SORT_NAME sortindirect
#define SORT_TYPE1 int
#define SORT_TYPE2 int
#include "nauty/sorttemplates.c"

#define SORT_OF_SORT 1
#define SORT_NAME sortints
#define SORT_TYPE1 int
#include "nauty/sorttemplates.c"

#define SORT_OF_SORT 2
#define SORT_NAME sortweights
#define SORT_TYPE1 int
#define SORT_TYPE2 sg_weight
#include "nauty/sorttemplates.c"
#include "invariant.h"

/*****************************************************************************
*                                                                            *
*  preparemarks1(N) and preparemarks2(N)                                     *
*  make vmarks array large enough to mark 0..N-1 and such that               *
*  the next RESETMARKS command will work correctly                           *
*                                                                            *
*****************************************************************************/

#if !MAXN
static void
preparemarks1(size_t nn, short** vmark1, size_t* vmark1_sz, short* vmark1_val)
{
    size_t oldsize;
    short *oldpos;

    oldsize = *vmark1_sz;
    oldpos = *vmark1;
    DYNALLOC1(short,*vmark1,*vmark1_sz,nn,"preparemarks");
    if (*vmark1_sz != oldsize || *vmark1 != oldpos) *vmark1_val = 32000;
}
#endif

#if !MAXN
static void
preparemarks2(size_t nn, short** vmark2, size_t* vmark2_sz, short* vmark2_val)
{
    size_t oldsize;
    short *oldpos;

    oldsize = *vmark2_sz;
    oldpos  = *vmark2;
    DYNALLOC1(short,*vmark2,*vmark2_sz,nn,"preparemarks");
    if (*vmark2_sz != oldsize || *vmark2 != oldpos) *vmark2_val = 32000;
}
#endif


/*****************************************************************************
*                                                                            *
*  refine_sg(g,lab,ptn,level,numcells,count,active,code,m,n) performs a      *
*  refinement operation on the partition at the specified level of the       *
*  partition nest (lab,ptn).  *numcells is assumed to contain the number of  *
*  cells on input, and is updated.  The initial set of active cells (alpha   *
*  in the paper) is specified in the set active.  Precisely, x is in active  *
*  iff the cell starting at index x in lab is active.                        *
*  The resulting partition is equitable if active is correct (see the paper  *
*  and the Guide).                                                           *
*  *code is set to a value which depends on the fine detail of the           *
*  algorithm, but which is independent of the labelling of the graph.        *
*  count is used for work space.                                             *
*                                                                            *
*****************************************************************************/

bool
dynamic_refine_sg(graph *g, int *lab, int *ptn, int level, int *numcells,
          int *count, set *active, int *code, int m, int n, invariant* I, workspace* W)
{
    int i,j,k,l,v1,v2,v3,isplit;
    int w1,w2,w3;
    long longcode;
    int *d,*e;
    int size,bigsize,bigpos;
    int nactive,hitcells;
    int lj,di,splitv;
    boolean trivsplit;
    size_t *v,vi,ii;

    SG_VDE(g,v,d,e);

#if !MAXN
    DYNALLOC1(int,W->work1,W->work1_sz,n,"refine_sg");
    DYNALLOC1(int,W->work2,W->work2_sz,n,"refine_sg");
    DYNALLOC1(int,W->work3,W->work3_sz,n,"refine_sg");
    DYNALLOC1(int,W->work4,W->work4_sz,n,"refine_sg");
#endif
#define CELLSTART W->work1
#define ACTIVE    W->work2
#define HITS      W->work3
#define HITCELL   W->work4

    PREPAREMARKS1(n, &W->vmark1, &W->vmark1_sz, &W->vmark1_val);
    PREPAREMARKS2(n, &W->vmark2, &W->vmark2_sz, &W->vmark2_val);

    longcode = *numcells;

    /* Set ACTIVE[0..nactive-1] = queue of active cell starts */

    nactive = 0;
    for (i = -1; (i = nextelement(active,m,i)) >= 0;)
        ACTIVE[nactive++] = i;

    if (nactive == 0)
    {
        *code = CLEANUP(longcode);
        return true;
    }

    /* Set CELLSTART[i] = starting point in lab[] of nontrivial cell
   containing i, or n if i is a singleton */

    for (i = 0; i < n; )
    {
        /* Just here, i is a cell starting position */
        if (ptn[i] <= level)
        {
            CELLSTART[lab[i]] = n;
            ++i;
        }
        else
        {
            j = i;
            do
            {
                CELLSTART[lab[i]] = j;
            } while (ptn[i++] > level);
        }
    }

    if (level <= 2 && nactive == 1 && ptn[ACTIVE[0]] <= level
        && *numcells <= n/8)
    {
        isplit = ACTIVE[--nactive];
        DELELEMENT(active,isplit);

        distvals((sparsegraph*)g,lab[isplit],HITS,n);

        for (v1 = 0; v1 < n; )
        {
            if (ptn[v1] <= level)
            {
                ++v1;
                continue;
            }

            longcode = MASH(longcode,v1);
            if(!I->write_top_and_compare((v1))) {return false;}
            w1 = HITS[lab[v1]];

            v2 = v1+1;
            while (ptn[v2-1] > level && HITS[lab[v2]] == w1) ++v2;

            if (ptn[v2-1] <= level)
            {
                v1 = v2;
                continue;
            }

            w2 = NAUTY_INFINITY;
            v3 = j = v2;

            do
            {
                lj = lab[j];
                w3 = HITS[lj];
                if (w3 == w1)
                {
                    lab[j] = lab[v3];
                    lab[v3] = lab[v2];
                    lab[v2] = lj;
                    ++v2;
                    ++v3;
                }
                else if (w3 == w2)
                {
                    lab[j] = lab[v3];
                    lab[v3] = lj;
                    ++v3;
                }
                else if (w3 < w1)
                {
                    lab[j] = lab[v2];
                    lab[v2] = lab[v1];
                    lab[v1] = lj;
                    v3 = v2 + 1;
                    v2 = v1 + 1;
                    w2 = w1;
                    w1 = w3;
                }
                else if (w3 < w2)
                {
                    lab[j] = lab[v2];
                    lab[v2] = lj;
                    v3 = v2 + 1;
                    w2 = w3;
                }
            } while (ptn[j++] > level);

            longcode = MASH(longcode,w2);
            if(!I->write_top_and_compare((w2))) {return false;}
            longcode = MASH(longcode,v2);
            if(!I->write_top_and_compare((v2))) {return false;}
            if (j != v2)   /* At least two fragments
                                * v1..v2-1 = w1; v2..v3-1 = w2  */
            {
                if (v2 == v1+1)
                    CELLSTART[lab[v1]] = n;

                if (v3 == v2+1)
                    CELLSTART[lab[v2]] = n;
                else
                    for (k = v2; k < v3; ++k)
                        CELLSTART[lab[k]] = v2;
                ++*numcells;
                ptn[v2-1] = level;

                if (j == v3)
                {
                    /* Two fragments only */
                    if (v2-v1 <= v3-v2 && !ISELEMENT(active,v1))
                    {
                        ADDELEMENT(active,v1);
                        ACTIVE[nactive++] = v1;
                    }
                    else
                    {
                        ADDELEMENT(active,v2);
                        ACTIVE[nactive++] = v2;
                    }
                }
                else
                {
                    /* Extra fragments: v3..j-1 > w2 */
                    sortindirect(lab+v3,HITS,j-v3);
                    ACTIVE[nactive++] = v2;
                    ADDELEMENT(active,v2);
                    if (v2-v1 >= v3-v2)
                    {
                        bigpos = -1;
                        bigsize = v2-v1;
                    }
                    else
                    {
                        bigpos = nactive-1;
                        bigsize = v3-v2;
                    }
                    for (k = v3-1; k < j-1;)
                    {
                        ptn[k] = level;
                        longcode = MASH(longcode,k);
                        if(!I->write_top_and_compare((k))) {return false;}
                        ++*numcells;
                        l = k+1;
                        ADDELEMENT(active,l);
                        ACTIVE[nactive++] = l;
                        w3 = HITS[lab[l]];
                        for (k = l; k < j-1
                                    && HITS[lab[k+1]] == w3; ++k)
                            CELLSTART[lab[k+1]] = l;
                        size = k-l+1;
                        if (size == 1)
                            CELLSTART[lab[l]] = n;
                        else
                        {
                            CELLSTART[lab[l]] = l;
                            if (size > bigsize)
                            {
                                bigsize = size;
                                bigpos = nactive-1;
                            }
                        }
                    }

                    if (bigpos >= 0 && !ISELEMENT(active,v1))
                    {
                        longcode = MASH(longcode,bigpos);
                        if(!I->write_top_and_compare((bigpos))) {return false;}
                        DELELEMENT(active,ACTIVE[bigpos]);
                        ADDELEMENT(active,v1);
                        ACTIVE[bigpos] = v1;
                    }
                }
            }
            v1 = j;
        }
    }

    /* Iterate until complete */
    while (nactive > 0 && *numcells < n)
    {
        for (i = 0; i < nactive && i < 10; ++i)
            if (ptn[ACTIVE[i]] <= level) break;

        if (i < nactive && i < 10)
        {
            trivsplit = TRUE;
            isplit = ACTIVE[i];
            ACTIVE[i] = ACTIVE[--nactive];
        }
        else
        {
            isplit = ACTIVE[--nactive];
            trivsplit = ptn[isplit] <= level;
        }

        DELELEMENT(active,isplit);
        longcode = MASH(longcode,isplit);
        if(!I->write_top_and_compare((isplit))) {return false;}

        if (trivsplit)
        {
            RESETMARKS1;
            RESETMARKS2;
            hitcells = 0;
            splitv = lab[isplit];
            vi = v[splitv];
            di = d[splitv];
            for (ii = 0; ii < di; ++ii)
            {
                j = e[vi+ii];
                MARK2(j);
                k = CELLSTART[j];
                if (k != n && ISNOTMARKED1(k))
                {
                    MARK1(k);
                    HITCELL[hitcells++] = k;
                }
            }

            if (hitcells > 1) sortints(HITCELL,hitcells);
            longcode = MASH(longcode,hitcells);
            if(!I->write_top_and_compare((hitcells))) {return false;}

            /* divide cells according to which vertices are hit */

            for (i = 0; i < hitcells; ++i)
            {
                j = v1 = v2 = HITCELL[i];
                longcode = MASH(longcode,v2);
                if(!I->write_top_and_compare((v2))) {return false;}
                k = 0;
                do
                {
                    lj = lab[j];
                    if (ISMARKED2(lj))
                        HITS[k++] = lj;
                    else
                        lab[v2++] = lj;
                } while (ptn[j++] > level);

                longcode = MASH(longcode,k);
                if(!I->write_top_and_compare((k))) {return false;}
                v3 = v2;
                while (--k >= 0)
                {
                    j = HITS[k];
                    CELLSTART[j] = v2;
                    lab[v3++] = j;
                }

                if (v2 != v3 && v2 != v1)
                {
                    ++*numcells;
                    if (v2 == v1+1) CELLSTART[lab[v1]] = n;
                    if (v3 == v2+1) CELLSTART[lab[v2]] = n;
                    ptn[v2-1] = level;
                    longcode = MASH(longcode,v2);
                    if(!I->write_top_and_compare((v2))) {return false;}
                    if (v2-v1 <= v3-v2 && !ISELEMENT(active,v1))
                    {
                        ADDELEMENT(active,v1);
                        ACTIVE[nactive++] = v1;
                    }
                    else
                    {
                        ADDELEMENT(active,v2);
                        ACTIVE[nactive++] = v2;
                    }
                }
            }
        }
        else  /* non-trivial splitting */
        {
            /* isplit is the start of the splitting cell.
               Set HITS[i] = hits of i for i in non-trivial cells,
               HITCELL[0..hitcells-1] = starts of hit non-trivial cells */

            RESETMARKS1;
            hitcells = 0;
            do
            {
                vi = v[lab[isplit]];
                di = d[lab[isplit]];
                for (ii = 0; ii < di; ++ii)
                {
                    j = e[vi+ii];
                    k = CELLSTART[j];
                    if (k != n)
                    {
                        if (ISNOTMARKED1(k))
                        {
                            MARK1(k);
                            HITCELL[hitcells++] = k;
                            do
                            {
                                HITS[lab[k]] = 0;
                            } while (ptn[k++] > level);
                        }
                        ++HITS[j];
                    }
                }
            } while (ptn[isplit++] > level);

            if (hitcells > 1) sortints(HITCELL,hitcells);

            /* divide cells according to hit counts */

            longcode = MASH(longcode,hitcells);
            if(!I->write_top_and_compare((hitcells))) {return false;}
            for (i = 0; i < hitcells; ++i)
            {
                v1 = HITCELL[i];
                w1 = HITS[lab[v1]];
                longcode = MASH(longcode,v1);
                if(!I->write_top_and_compare((v1))) {return false;}

                v2 = v1+1;
                while (ptn[v2-1] > level && HITS[lab[v2]] == w1) ++v2;

                if (ptn[v2-1] <= level) continue;
                w2 = NAUTY_INFINITY;
                v3 = j = v2;

                do
                {
                    lj = lab[j];
                    w3 = HITS[lj];
                    if (w3 == w1)
                    {
                        lab[j] = lab[v3];
                        lab[v3] = lab[v2];
                        lab[v2] = lj;
                        ++v2;
                        ++v3;
                    }
                    else if (w3 == w2)
                    {
                        lab[j] = lab[v3];
                        lab[v3] = lj;
                        ++v3;
                    }
                    else if (w3 < w1)
                    {
                        lab[j] = lab[v2];
                        lab[v2] = lab[v1];
                        lab[v1] = lj;
                        v3 = v2 + 1;
                        v2 = v1 + 1;
                        w2 = w1;
                        w1 = w3;
                    }
                    else if (w3 < w2)
                    {
                        lab[j] = lab[v2];
                        lab[v2] = lj;
                        v3 = v2 + 1;
                        w2 = w3;
                    }
                } while (ptn[j++] > level);

                longcode = MASH(longcode,w1);
                if(!I->write_top_and_compare((w1))) {return false;}
                longcode = MASH(longcode,v2);
                if(!I->write_top_and_compare((v2))) {return false;}
                if (j != v2)   /* At least two fragments
                                * v1..v2-1 = w1; v2..v3-1 = w2  */
                {
                    if (v2 == v1+1)
                        CELLSTART[lab[v1]] = n;

                    if (v3 == v2+1)
                        CELLSTART[lab[v2]] = n;
                    else
                        for (k = v2; k < v3; ++k)
                            CELLSTART[lab[k]] = v2;
                    ++*numcells;
                    ptn[v2-1] = level;

                    if (j == v3)
                    {
                        /* Two fragments only */
                        if (v2-v1 <= v3-v2 && !ISELEMENT(active,v1))
                        {
                            ADDELEMENT(active,v1);
                            ACTIVE[nactive++] = v1;
                        }
                        else
                        {
                            ADDELEMENT(active,v2);
                            ACTIVE[nactive++] = v2;
                        }
                    }
                    else
                    {
                        /* Extra fragments: v3..j-1 > w2 */
                        longcode = MASH(longcode,v3);
                        if(!I->write_top_and_compare((v3))) {return false;}
                        sortindirect(lab+v3,HITS,j-v3);
                        ACTIVE[nactive++] = v2;
                        ADDELEMENT(active,v2);
                        if (v2-v1 >= v3-v2)
                        {
                            bigpos = -1;
                            bigsize = v2-v1;
                        }
                        else
                        {
                            bigpos = nactive-1;
                            bigsize = v3-v2;
                            longcode = MASH(longcode,bigsize);
                            if(!I->write_top_and_compare((bigsize))) {return false;}
                        }
                        for (k = v3-1; k < j-1;)
                        {
                            ptn[k] = level;
                            ++*numcells;
                            l = k+1;
                            ADDELEMENT(active,l);
                            ACTIVE[nactive++] = l;
                            w3 = HITS[lab[l]];
                            longcode = MASH(longcode,w3);
                            if(!I->write_top_and_compare((w3))) {return false;}
                            for (k = l; k < j-1
                                        && HITS[lab[k+1]] == w3; ++k)
                                CELLSTART[lab[k+1]] = l;
                            size = k-l+1;
                            if (size == 1)
                                CELLSTART[lab[l]] = n;
                            else
                            {
                                CELLSTART[lab[l]] = l;
                                if (size > bigsize)
                                {
                                    bigsize = size;
                                    bigpos = nactive-1;
                                }
                            }
                        }

                        if (bigpos >= 0 && !ISELEMENT(active,v1))
                        {
                            DELELEMENT(active,ACTIVE[bigpos]);
                            ADDELEMENT(active,v1);
                            ACTIVE[bigpos] = v1;
                        }
                    }
                }
            }
        }
    }

    longcode = MASH(longcode,*numcells);
    if(!I->write_top_and_compare((*numcells))) {return false;}
    *code = CLEANUP(longcode);
    return true;
}