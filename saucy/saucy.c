/*
 * saucy.c
 * Searching for Automorphisms in Underlying CNF, yes?
 *
 * by Paul T. Darga <pdarga@umich.edu>
 * and Mark Liffiton <liffiton@umich.edu>
 * and Hadi Katebi <hadik@eecs.umich.edu>
 * 
 * Copyright (C) 2004, The Regents of the University of Michigan
 * See the LICENSE file for details.
 */

#include <stdlib.h> /* malloc, calloc, and free */
#include <string.h> /* memcpy */

#include "saucy.h"

struct coloring {
	int *lab;        /* Labelling of objects */
	int *unlab;      /* Inverse of lab */
	int *cfront;     /* Pointer to front of cells */
	int *clen;       /* Length of cells (defined for cfront's) */
};

struct saucy {
	/* Graph data */
	int n;           /* Size of domain */
	const int *adj;  /* Neighbors of k: edg[adj[k]]..edg[adj[k+1]] */
	const int *edg;  /* Actual neighbor data */
	const int *dadj; /* Fanin neighbor indices, for digraphs */
	const int *dedg; /* Fanin neighbor data, for digraphs */
	void *arg;       /* Opaque client data */

	/* Coloring data */
	struct coloring left, right;
	int *nextnon;    /* Forward next-nonsingleton pointers */ 
	int *prevnon;    /* Backward next-nonsingleton pointers */

	/* Refinement: inducers */
	char *indmark;   /* Induce marks */
	int *ninduce;    /* Nonsingletons that might induce refinement */
	int *sinduce;    /* Singletons that might induce refinement */
	int nninduce;    /* Size of ninduce stack */
	int nsinduce;    /* Size of sinduce stack */

	/* Refinement: marked cells */
	int *clist;      /* List of cells marked for refining */
	int csize;       /* Number of cells in clist */

	/* Refinement: workspace */
	char *stuff;     /* Bit vector, but one char per bit */
	int *ccount;     /* Number of connections to refining cell */
	int *bucket;     /* Workspace */
	int *count;      /* Num vertices with same adj count to ref cell */
	int *junk;       /* More workspace */
	int *gamma;      /* Working permutation */
	int *conncnts;   /* Connection counts for cell fronts */

	/* Search data */
	int lev;         /* Current search tree level */
	int anc;         /* Level of greatest common ancestor with zeta */
	int *anctar;     /* Copy of target cell at anc */
	int kanctar;     /* Location within anctar to iterate from */
	int *start;      /* Location of target at each level */
	int indmin;      /* Used for group size computation */
	int match;       /* Have we not diverged from previous left? */

	/* Search: orbit partition */
	int *theta;      /* Running approximation of orbit partition */
	int *thsize;     /* Size of cells in theta, defined for mcrs */
	int *thnext;     /* Next rep in list (circular list) */
	int *thprev;     /* Previous rep in list */
	int *threp;      /* First rep for a given cell front */
	int *thfront;    /* The cell front associated with this rep */

	/* Search: split record */
	int *splitwho;   /* List of where splits occurred */
	int *splitfrom;  /* List of cells which were split */
	int *splitlev;   /* Where splitwho/from begins for each level */
	int nsplits;     /* Number of splits at this point */

	/* Search: differences from leftmost */
	char *diffmark;  /* Marked for diff labels */
	int *diffs;      /* List of diff labels */
	int *difflev;    /* How many labels diffed at each level */
	int ndiffs;      /* Current number of diffs */
	int *undifflev;  /* How many diff labels fixed at each level */
	int nundiffs;    /* Current number of diffs in singletons (fixed) */
	int *unsupp;     /* Inverted diff array */
	int *specmin;    /* Speculated mappings */
	int *pairs;      /* Not-undiffed diffs that can make two-cycles */
	int *unpairs;    /* Indices into pairs */
	int npairs;      /* Number of such pairs */
	int *diffnons;   /* Diffs that haven't been undiffed */
	int *undiffnons; /* Inverse of that */
	int ndiffnons;   /* Number of such diffs */

	/* Polymorphic functions */
	saucy_consumer *consumer;
	int (*split)(struct saucy *, struct coloring *, int, int);
	int (*is_automorphism)(struct saucy *);
	int (*ref_singleton)(struct saucy *, struct coloring *, int);
	int (*ref_nonsingle)(struct saucy *, struct coloring *, int);

	 /* Statistics structure */
	struct saucy_stats *stats;
};

static int
array_find_min(const int *a, int n)
{
	const int *start = a, *end = a + n, *min = a;
	while (++a != end) {
		if (*a < *min) min = a;
	}
	return min - start;
}

static void
swap(int *a, int x, int y)
{
	int tmp = a[x];
	a[x] = a[y];
	a[y] = tmp;
}

static void
sift_up(int *a, int k)
{
	int p;
	do {
		p = k / 2;
		if (a[k] <= a[p]) {
			return;
		}
		else {
			swap(a, k, p);
			k = p;
		}
	} while (k > 1);
}

static void
sift_down(int *a, int n)
{
	int p = 1, k = 2;
	while (k <= n) {
		if (k < n && a[k] < a[k+1]) ++k;
		if (a[p] < a[k]) {
			swap(a, p, k);
			p = k;
			k = 2 * p;
		}
		else {
			return;
		}
	}
}

static void
heap_sort(int *a, int n)
{
	int i;
	for (i = 1; i < n; ++i) {
		sift_up(a-1, i+1);
	}
	--i;
	while (i > 0) {
		swap(a, 0, i);
		sift_down(a-1, i--);
	}
}

static void
insertion_sort(int *a, int n)
{
	int i, j, k;
	for (i = 1; i < n; ++i) {
		k = a[i];
		for (j = i; j > 0 && a[j-1] > k; --j) {
			a[j] = a[j-1];
		}
		a[j] = k;
	}
}

static int
partition(int *a, int n, int m)
{
	int f = 0, b = n;
	for (;;) {
		while (a[f] <= m) ++f;
		do  --b; while (m <= a[b]);
		if (f < b) {
			swap(a, f, b);
			++f;
		}
		else break;
	}
	return f;
}

static int
log_base2(int n)
{
	int k = 0;
	while (n > 1) {
		++k;
		n >>= 1;
	}
	return k;
}

static int
median(int a, int b, int c)
{
	if (a <= b) {
		if (b <= c) return b;
		if (a <= c) return c;
		return a;
	}
	else {
		if (a <= c) return a;
		if (b <= c) return c;
		return b;
	}
}

static void
introsort_loop(int *a, int n, int lim)
{
	int p;
	while (n > 16) {
		if (lim == 0) {
			heap_sort(a, n);
			return;
		}
		--lim;
		p = partition(a, n, median(a[0], a[n/2], a[n-1]));
		introsort_loop(a + p, n - p, lim);
		n = p;
	}
}

static void
introsort(int *a, int n)
{
	introsort_loop(a, n, 2 * log_base2(n));
	insertion_sort(a, n);
}

static int
do_find_min(struct coloring *c, int t)
{
	return array_find_min(c->lab + t, c->clen[t] + 1) + t;
}

static int
find_min(struct saucy *s, int t)
{
	return do_find_min(&s->right, t);
}

static void
set_label(struct coloring *c, int index, int value)
{
	c->lab[index] = value;
	c->unlab[value] = index;
}

static void
swap_labels(struct coloring *c, int a, int b)
{
	int tmp = c->lab[a];
	set_label(c, a, c->lab[b]);
	set_label(c, b, tmp);
}

static void
move_to_back(struct saucy *s, struct coloring *c, int k)
{
	int cf = c->cfront[k];
	int cb = cf + c->clen[cf];
	int offset = s->conncnts[cf]++;

	/* Move this connected label to the back of its cell */
	swap_labels(c, cb - offset, c->unlab[k]);

	/* Add it to the cell list if it's the first one swapped */
	if (!offset) s->clist[s->csize++] = cf;
}

static void
data_mark(struct saucy *s, struct coloring *c, int k)
{
	int cf = c->cfront[k];

	/* Move connects to the back of nonsingletons */
	if (c->clen[cf]) move_to_back(s, c, k);
}

static void
data_count(struct saucy *s, struct coloring *c, int k)
{
	int cf = c->cfront[k];

	/* Move to back and count the number of connections */
	if (c->clen[cf] && !s->ccount[k]++) move_to_back(s, c, k);
}

static int
check_mapping(struct saucy *s, const int *adj, const int *edg, int k)
{
	int i, gk, ret = 1;

	/* Mark gamma of neighbors */
	for (i = adj[k]; i != adj[k+1]; ++i) {
		s->stuff[s->gamma[edg[i]]] = 1;
	}

	/* Check neighbors of gamma */
	gk = s->gamma[k];
	for (i = adj[gk]; ret && i != adj[gk+1]; ++i) {
		ret = s->stuff[edg[i]];
	}

	/* Clear out bit vector before we leave */
	for (i = adj[k]; i != adj[k+1]; ++i) {
		s->stuff[s->gamma[edg[i]]] = 0;
	}

	return ret;
}

static int
is_undirected_automorphism(struct saucy *s)
{
	int i, j;

	for (i = 0; i < s->ndiffs; ++i) {
		j = s->unsupp[i];
		if (!check_mapping(s, s->adj, s->edg, j)) return 0;
	}
	return 1;
}

static int
is_directed_automorphism(struct saucy *s)
{
	int i, j;

	for (i = 0; i < s->ndiffs; ++i) {
		j = s->unsupp[i];
		if (!check_mapping(s, s->adj, s->edg, j)) return 0;
		if (!check_mapping(s, s->dadj, s->dedg, j)) return 0;
	}
	return 1;
}

static void
add_induce(struct saucy *s, struct coloring *c, int who)
{
	if (!c->clen[who]) {
		s->sinduce[s->nsinduce++] = who;
	}
	else {
		s->ninduce[s->nninduce++] = who;
	}
	s->indmark[who] = 1;
}

static void
fix_fronts(struct coloring *c, int cf, int ff)
{
	int i, end = cf + c->clen[cf];
	for (i = ff; i <= end; ++i) {
		c->cfront[c->lab[i]] = cf;
	}
}

static void
array_indirect_sort(int *a, const int *b, int n)
{
	int h, i, j, k;

	/* Shell sort, as implemented in nauty, (C) Brendan McKay */
	j = n / 3;
	h = 1;
	do { h = 3 * h + 1; } while (h < j);

	do {
		for (i = h; i < n; ++i) {
			k = a[i];
			for (j = i; b[a[j-h]] > b[k]; ) {
				a[j] = a[j-h];
				if ((j -= h) < h) break;
			}
			a[j] = k;
		}
		h /= 3;
	} while (h > 0);
}

static int
at_terminal(struct saucy *s)
{
	return s->nsplits == s->n;
}

static void
add_diffnon(struct saucy *s, int k)
{
	/* Only add if we're in a consistent state */
	if (s->ndiffnons == -1) return;

	s->undiffnons[k] = s->ndiffnons;
	s->diffnons[s->ndiffnons++] = k;
}

static void
remove_diffnon(struct saucy *s, int k)
{
	int j;

	if (s->undiffnons[k] == -1) return;

	j = s->diffnons[--s->ndiffnons];
	s->diffnons[s->undiffnons[k]] = j;
	s->undiffnons[j] = s->undiffnons[k];

	s->undiffnons[k] = -1;
}

static void
add_diff(struct saucy *s, int k)
{
	if (!s->diffmark[k]) {
		s->diffmark[k] = 1;
		s->diffs[s->ndiffs++] = k;
		add_diffnon(s, k);
	}
}

static int
is_a_pair(struct saucy *s, int k)
{
	return s->unpairs[k] != -1;
}

static int
in_cell_range(struct coloring *c, int ff, int cf)
{
	int cb = cf + c->clen[cf];
	return cf <= ff && ff <= cb;
}

static void
add_pair(struct saucy *s, int k)
{
	if (s->npairs != -1) {
		s->unpairs[k] = s->npairs;
		s->pairs[s->npairs++] = k;
	}
}

static void
eat_pair(struct saucy *s, int k)
{
	int j;
	j = s->pairs[--s->npairs];
	s->pairs[s->unpairs[k]] = j;
	s->unpairs[j] = s->unpairs[k];
	s->unpairs[k] = -1;
}

static void
pick_all_the_pairs(struct saucy *s)
{
	int i;
	for (i = 0; i < s->npairs; ++i) {
		s->unpairs[s->pairs[i]] = -1;
	}
	s->npairs = 0;
}

static void
clear_undiffnons(struct saucy *s)
{
	int i;
	for (i = 0 ; i < s->ndiffnons ; ++i) {
		s->undiffnons[s->diffnons[i]] = -1;
	}
}

static void
fix_diff_singleton(struct saucy *s, int cf)
{
	int r = s->right.lab[cf];
	int l = s->left.lab[cf];
	int rcfl;

	if (!s->right.clen[cf] && r != l) {

		/* Make sure diff is marked */
		add_diff(s, r);

		/* It is now undiffed since it is singleton */
		++s->nundiffs;
		remove_diffnon(s, r);

		/* Mark the other if not singleton already */
		rcfl = s->right.cfront[l];
		if (s->right.clen[rcfl]) {
			add_diff(s, l);

			/* Check for pairs */
			if (in_cell_range(&s->right, s->left.unlab[r], rcfl)) {
				add_pair(s, l);
			}
		}
		/* Otherwise we might be eating a pair */
		else if (is_a_pair(s, r)) {
			eat_pair(s, r);
		}
	}
}

static void
fix_diff_subtract(struct saucy *s, int cf, const int *a, const int *b)
{
	int i, k;
	int cb = cf + s->right.clen[cf];

	/* Mark the contents of the first set */
	for (i = cf; i <= cb; ++i) {
		s->stuff[a[i]] = 1;
	}

	/* Add elements from second set not present in the first */
	for (i = cf; i <= cb; ++i) {
		k = b[i];
		if (!s->stuff[k]) add_diff(s, k);
	}

	/* Clear the marks of the first set */
	for (i = cf; i <= cb; ++i) {
		s->stuff[a[i]] = 0;
	}
}

static void
fix_diffs(struct saucy *s, int cf, int ff)
{
	int min;

	/* Check for singleton cases in both cells */
	fix_diff_singleton(s, cf);
	fix_diff_singleton(s, ff);

	/* If they're both nonsingleton, do subtraction on smaller */
	if (s->right.clen[cf] && s->right.clen[ff]) {
		min = s->right.clen[cf] < s->right.clen[ff] ? cf : ff;
		fix_diff_subtract(s, min, s->left.lab, s->right.lab);
		fix_diff_subtract(s, min, s->right.lab, s->left.lab);
	}
}

static void
split_color(struct coloring *c, int cf, int ff)
{
	int cb, fb;

	/* Fix lengths */
	fb = ff - 1;
	cb = cf + c->clen[cf];
	c->clen[cf] = fb - cf;
	c->clen[ff] = cb - ff;

	/* Fix cell front pointers */
	fix_fronts(c, ff, ff);
}

static void
split_common(struct saucy *s, struct coloring *c, int cf, int ff)
{
	split_color(c, cf, ff);

	/* Add to refinement */
	if (s->indmark[cf] || c->clen[ff] < c->clen[cf]) {
		add_induce(s, c, ff);
	}
	else {
		add_induce(s, c, cf);
	}
}

static int
split_left(struct saucy *s, struct coloring *c, int cf, int ff)
{
	/* Record the split */
	s->splitwho[s->nsplits] = ff;
	s->splitfrom[s->nsplits] = cf;
	++s->nsplits;

	/* Do common splitting tasks */
	split_common(s, c, cf, ff);

	/* Always succeeds */
	return 1;
}

static int
split_init(struct saucy *s, struct coloring *c, int cf, int ff)
{
	split_left(s, c, cf, ff);

	/* Maintain nonsingleton list for finding new targets */
	if (c->clen[ff]) {
		s->prevnon[s->nextnon[cf]] = ff;
		s->nextnon[ff] = s->nextnon[cf];
		s->prevnon[ff] = cf;
		s->nextnon[cf] = ff;
	}
	if (!c->clen[cf]) {
		s->nextnon[s->prevnon[cf]] = s->nextnon[cf];
		s->prevnon[s->nextnon[cf]] = s->prevnon[cf];
	}

	/* Always succeeds */
	return 1;
}

static int
split_other(struct saucy *s, struct coloring *c, int cf, int ff)
{
	int k = s->nsplits;

	/* Verify the split with init */
	if (s->splitwho[k] != ff || s->splitfrom[k] != cf
			|| k >= s->splitlev[s->lev]) {
		return 0;
	}
	++s->nsplits;

	/* Do common splitting tasks */
	split_common(s, c, cf, ff);

	/* Fix differences with init */
	fix_diffs(s, cf, ff);

	/* If we got this far we succeeded */
	return 1;
}

static int
refine_cell(struct saucy *s, struct coloring *c,
	int (*refine)(struct saucy *, struct coloring *, int))
{
	int i, cf, ret = 1;

	/*
	 * The connected list must be consistent.  This is for
	 * detecting mappings across nodes at a given level.  However,
	 * at the root of the tree, we never have to map with another
	 * node, so we lack this consistency constraint in that case.
	 */
	if (s->lev > 1) introsort(s->clist, s->csize);

	/* Now iterate over the marked cells */
	for (i = 0; ret && i < s->csize; ++i) {
		cf = s->clist[i];
		ret = refine(s, c, cf);
	}

	/* Clear the connected marks */
	for (i = 0; i < s->csize; ++i) {
		cf = s->clist[i];
		s->conncnts[cf] = 0;
	}
	s->csize = 0;
	return ret;
}

static int
maybe_split(struct saucy *s, struct coloring *c, int cf, int ff)
{
	return cf == ff ? 1 : s->split(s, c, cf, ff);
}

static int
ref_single_cell(struct saucy *s, struct coloring *c, int cf)
{
	int zcnt = c->clen[cf] + 1 - s->conncnts[cf];
	return maybe_split(s, c, cf, cf + zcnt);
}

static int
ref_singleton(struct saucy *s, struct coloring *c,
	const int *adj, const int *edg, int cf)
{
	int i, k = c->lab[cf];

	/* Find the cells we're connected to, and mark our neighbors */
	for (i = adj[k]; i != adj[k+1]; ++i) {
		data_mark(s, c, edg[i]);
	}

	/* Refine the cells we're connected to */
	return refine_cell(s, c, ref_single_cell);
}

static int
ref_singleton_directed(struct saucy *s, struct coloring *c, int cf)
{
	return ref_singleton(s, c, s->adj, s->edg, cf)
		&& ref_singleton(s, c, s->dadj, s->dedg, cf);
}

static int
ref_singleton_undirected(struct saucy *s, struct coloring *c, int cf)
{
	return ref_singleton(s, c, s->adj, s->edg, cf);
}

static int
ref_nonsingle_cell(struct saucy *s, struct coloring *c, int cf)
{
	int cnt, i, cb, nzf, ff, fb, bmin, bmax;

	/* Find the front and back */
	cb = cf + c->clen[cf];
	nzf = cb - s->conncnts[cf] + 1;

	/* Prepare the buckets */
	ff = nzf;
	cnt = s->ccount[c->lab[ff]];
	s->count[ff] = bmin = bmax = cnt;
	s->bucket[cnt] = 1;

	/* Iterate through the rest of the vertices */
	while (++ff <= cb) {
		cnt = s->ccount[c->lab[ff]];

		/* Initialize intermediate buckets */
		while (bmin > cnt) s->bucket[--bmin] = 0;
		while (bmax < cnt) s->bucket[++bmax] = 0;

		/* Mark this count */
		++s->bucket[cnt];
		s->count[ff] = cnt;
	}

	/* If they all had the same count, bail */
	if (bmin == bmax && cf == nzf) return 1;
	ff = fb = nzf;

	/* Calculate bucket locations, sizes */
	for (i = bmin; i <= bmax; ++i, ff = fb) {
		if (!s->bucket[i]) continue;
		fb = ff + s->bucket[i];
		s->bucket[i] = fb;
	}

	/* Repair the partition nest */
	for (i = nzf; i <= cb; ++i) {
		s->junk[--s->bucket[s->count[i]]] = c->lab[i];
	}
	for (i = nzf; i <= cb; ++i) {
		set_label(c, i, s->junk[i]);
	}

	/* Split; induce */
	for (i = bmax; i > bmin; --i) {
		ff = s->bucket[i];
		if (ff && !s->split(s, c, cf, ff)) return 0;
	}

	/* If there was a zero area, then there's one more cell */
	return maybe_split(s, c, cf, s->bucket[bmin]);
}

static int
ref_nonsingle(struct saucy *s, struct coloring *c,
	const int *adj, const int *edg, int cf)
{
	int i, j, k, ret;
	const int cb = cf + c->clen[cf];
	const int size = cb - cf + 1;

	/* Double check for nonsingles which became singles later */
	if (cf == cb) {
		return ref_singleton(s, c, adj, edg, cf);
	}

	/* Establish connected list */
	memcpy(s->junk, c->lab + cf, size * sizeof(int));
	for (i = 0; i < size; ++i) {
		k = s->junk[i];
		for (j = adj[k]; j != adj[k+1]; ++j) {
			data_count(s, c, edg[j]);
		}
	}

	/* Refine the cells we're connected to */
	ret = refine_cell(s, c, ref_nonsingle_cell);

	/* Clear the counts; use lab because junk was overwritten */
	for (i = cf; i <= cb; ++i) {
		k = c->lab[i];
		for (j = adj[k]; j != adj[k+1]; ++j) {
			s->ccount[edg[j]] = 0;
		}
	}

	return ret;
}

static int
ref_nonsingle_directed(struct saucy *s, struct coloring *c, int cf)
{
	return ref_nonsingle(s, c, s->adj, s->edg, cf)
		&& ref_nonsingle(s, c, s->dadj, s->dedg, cf);
}

static int
ref_nonsingle_undirected(struct saucy *s, struct coloring *c, int cf)
{
	return ref_nonsingle(s, c, s->adj, s->edg, cf);
}

static void
clear_refine(struct saucy *s)
{
	int i;
	for (i = 0; i < s->nninduce; ++i) {
		s->indmark[s->ninduce[i]] = 0;
	}
	for (i = 0; i < s->nsinduce; ++i) {
		s->indmark[s->sinduce[i]] = 0;
	}
	s->nninduce = s->nsinduce = 0;
}

static int
refine(struct saucy *s, struct coloring *c)
{
	int front;

	/* Keep going until refinement stops */
	while (1) {

		/* If discrete, bail */
		if (at_terminal(s)) {
			clear_refine(s);
			return 1;
		};

		/* Look for something else to refine on */
		if (s->nsinduce) {
			front = s->sinduce[--s->nsinduce];
			s->indmark[front] = 0;
			if (!s->ref_singleton(s, c, front)) break;
		}
		else if (s->nninduce) {
			front = s->ninduce[--s->nninduce];
			s->indmark[front] = 0;
			if (!s->ref_nonsingle(s, c, front)) break;
		}
		else {
			return 1;
		};
	}

	clear_refine(s);
	return 0;
}

static int
descend(struct saucy *s, struct coloring *c, int target, int min)
{
	int back = target + c->clen[target];
	int ret;

	/* Count this node */
	++s->stats->nodes;

	/* Move the minimum label to the back */
	swap_labels(c, min, back);

	/* Split the cell */
	s->difflev[s->lev] = s->ndiffs;
	s->undifflev[s->lev] = s->nundiffs;
	++s->lev;
	s->split(s, c, target, back);

	/* Now go and do some work */
	ret = refine(s, c);

	/* This is the new enhancement in saucy 3.0 */
	if (c == &s->right && ret) {
        	int i, j, v, sum1, sum2, xor1, xor2;
		for (i = s->nsplits - 1; i > s->splitlev[s->lev-1]; --i) {
			v = c->lab[s->splitwho[i]];
			sum1 = xor1 = 0;
			for (j = s->adj[v]; j < s->adj[v+1]; j++) {
				sum1 += c->cfront[s->edg[j]];
				xor1 ^= c->cfront[s->edg[j]];
			}
			v = s->left.lab[s->splitwho[i]];
			sum2 = xor2 = 0;
			for (j = s->adj[v]; j < s->adj[v+1]; j++) {
				sum2 += s->left.cfront[s->edg[j]];
				xor2 ^= s->left.cfront[s->edg[j]];
			}
			if ((sum1 != sum2) || (xor1 != xor2)) {
				ret = 0;
				break;
			}
			v = c->lab[s->splitfrom[i]];
			sum1 = xor1 = 0;
			for (j = s->adj[v]; j < s->adj[v+1]; j++) {
				sum1 += c->cfront[s->edg[j]];
				xor1 ^= c->cfront[s->edg[j]];
			}
			v = s->left.lab[s->splitfrom[i]];
			sum2 = xor2 = 0;
			for (j = s->adj[v]; j < s->adj[v+1]; j++) {
				sum2 += s->left.cfront[s->edg[j]];
				xor2 ^= s->left.cfront[s->edg[j]];
			}
			if ((sum1 != sum2) || (xor1 != xor2)) {
				ret = 0;
				break;
			}
		}
	}

	return ret;
}

static int
descend_leftmost(struct saucy *s)
{
	int target, min;

	/* Keep going until we're discrete */
	while (!at_terminal(s)) {
		target = s->nextnon[-1];
		min = target;
		s->start[s->lev] = target;
		s->splitlev[s->lev] = s->nsplits;
		if (!descend(s, &s->left, target, min)) return 0;
	}
	s->splitlev[s->lev] = s->n;
	return 1;
}

/*
 * If the remaining nonsingletons in this partition match exactly
 * those nonsingletons from the leftmost branch of the search tree,
 * then there is no point in continuing descent.
 */

static int
zeta_fixed(struct saucy *s)
{
	return s->ndiffs == s->nundiffs;
}

static void
select_decomposition(struct saucy *s, int *target, int *lmin, int *rmin)
{
	/* Both clens are equal; this clarifies the code a bit */
	const int *clen = s->left.clen;
	int i, cf, k;

	/*
	 * If there's a pair, use it.  pairs[0] should always work,
	 * but we use a checked loop instead because I'm not 100% sure
	 * I'm "unpairing" at every point I should be.
	 */
	for (i = 0; i < s->npairs; ++i) {
		k = s->pairs[i];
		*target = s->right.cfront[k];
		*lmin = s->left.unlab[s->right.lab[s->left.unlab[k]]];
		*rmin = s->right.unlab[k];

		if (clen[*target]
				&& in_cell_range(&s->left, *lmin, *target)
				&& in_cell_range(&s->right, *rmin, *target))
			return;
	}

	/* Diffnons is only consistent when there are no baddies */
	if (s->ndiffnons != -1) {
		*target = *lmin = *rmin = s->right.cfront[s->diffnons[0]];
		return;
	}

	/* Pick any old target cell and element */
	for (i = 0; i < s->ndiffs; ++i) {
		cf = s->right.cfront[s->diffs[i]];
		if (clen[cf]) {
			*lmin = *rmin = *target = cf;
			return;
		}
	}

	/* Shouldn't get here */
	abort();
}

static int
descend_left(struct saucy *s)
{
	int target, lmin, rmin;

	/* Check that we ended at the right spot */
	if (s->nsplits != s->splitlev[s->lev]) return 0;

	/* Keep going until we're discrete */
	while (!at_terminal(s) && !zeta_fixed(s)) {

		/* We can pick any old target cell and element */
		select_decomposition(s, &target, &lmin, &rmin);

		/* Check if we need to refine on the left */
		s->match = 0;
		s->start[s->lev] = target;
		s->split = split_left;
		descend(s, &s->left, target, lmin);
		s->splitlev[s->lev] = s->nsplits;
		s->split = split_other;
		--s->lev;
		s->nsplits = s->splitlev[s->lev];

		/* Now refine on the right and ensure matching */
		s->specmin[s->lev] = s->right.lab[rmin];
		if (!descend(s, &s->right, target, rmin)) return 0;
		if (s->nsplits != s->splitlev[s->lev]) return 0;
	}
	return 1;
}

static int
find_representative(int k, int *theta)
{
	int rep, tmp;

	/* Find the minimum cell representative */
	for (rep = k; rep != theta[rep]; rep = theta[rep]);

	/* Update all intermediaries */
	while (theta[k] != rep) {
		tmp = theta[k]; theta[k] = rep; k = tmp;
	}
	return rep;
}

static void
update_theta(struct saucy *s)
{
	int i, k, x, y, tmp;

	for (i = 0; i < s->ndiffs; ++i) {
		k = s->unsupp[i];
		x = find_representative(k, s->theta);
		y = find_representative(s->gamma[k], s->theta);

		if (x != y) {
			if (x > y) {
				tmp = x;
				x = y;
				y = tmp;
			}
			s->theta[y] = x;
			s->thsize[x] += s->thsize[y];

			s->thnext[s->thprev[y]] = s->thnext[y];
			s->thprev[s->thnext[y]] = s->thprev[y];
			s->threp[s->thfront[y]] = s->thnext[y];
		}
	}
}

static int
theta_prune(struct saucy *s)
{
	int start = s->start[s->lev];
	int label, rep, irep;

	irep = find_representative(s->indmin, s->theta);
	while (s->kanctar) {
		label = s->anctar[--s->kanctar];
		rep = find_representative(label, s->theta);
		if (rep == label && rep != irep) {
			return s->right.unlab[label] - start;
		}
	}
	return -1;
}

static int
orbit_prune(struct saucy *s)
{
	int i, label, fixed, min = -1;
	int k = s->start[s->lev];
	int size = s->right.clen[k] + 1;
	int *cell = s->right.lab + k;

	/* The previously fixed value */
	fixed = cell[size-1];

	/* Look for the next minimum cell representative */
	for (i = 0; i < size-1; ++i) {
		label = cell[i];

		/* Skip things we've already considered */
		if (label <= fixed) continue;

		/* Skip things that we'll consider later */
		if (min != -1 && label > cell[min]) continue;

		/* New candidate minimum */
		min = i;
	}

	return min;
}

static void
note_anctar_reps(struct saucy *s)
{
	int i, j, k, m, f, rep, tmp;

	/*
	 * Undo the previous level's splits along leftmost so that we
	 * join the appropriate lists of theta reps.
	 */
	for (i = s->splitlev[s->anc+1]-1; i >= s->splitlev[s->anc]; --i) {
		f = s->splitfrom[i];
		j = s->threp[f];
		k = s->threp[s->splitwho[i]];

		s->thnext[s->thprev[j]] = k;
		s->thnext[s->thprev[k]] = j;

		tmp = s->thprev[j];
		s->thprev[j] = s->thprev[k];
		s->thprev[k] = tmp;

		for (m = k; m != j; m = s->thnext[m]) {
			s->thfront[m] = f;
		}
	}

	/*
	 * Just copy over the target's reps and sort by cell size, in
	 * the hopes of trimming some otherwise redundant generators.
	 */
	s->kanctar = 0;
	s->anctar[s->kanctar++] = rep = s->threp[s->start[s->lev]];
	for (k = s->thnext[rep]; k != rep; k = s->thnext[k]) {
		s->anctar[s->kanctar++] = k;
	}
	array_indirect_sort(s->anctar, s->thsize, s->kanctar);
}

static void
multiply_index(struct saucy *s, int k)
{
	if ((s->stats->grpsize_base *= k) > 1e10) {
		s->stats->grpsize_base /= 1e10;
		s->stats->grpsize_exp += 10;
	}
}

static int
backtrack_leftmost(struct saucy *s)
{
	int rep = find_representative(s->indmin, s->theta);
	int repsize = s->thsize[rep];
	int min = -1;

	pick_all_the_pairs(s);
	clear_undiffnons(s);
	s->ndiffs = s->nundiffs = s->npairs = s->ndiffnons = 0;

	if (repsize != s->right.clen[s->start[s->lev]]+1) {
		min = theta_prune(s);
	}

	if (min == -1) {
		multiply_index(s, repsize);
	}

	return min;
}

static int
backtrack_other(struct saucy *s)
{
	int cf = s->start[s->lev];
	int cb = cf + s->right.clen[cf];
	int spec = s->specmin[s->lev];
	int min;

	/* Avoid using pairs until we get back to leftmost. */
	pick_all_the_pairs(s);

	clear_undiffnons(s);

	s->npairs = s->ndiffnons = -1;

	if (s->right.lab[cb] == spec) {
		min = find_min(s, cf);
		if (min == cb) {
			min = orbit_prune(s);
		}
		else {
			min -= cf;
		}
	}
	else {
		min = orbit_prune(s);
		if (min != -1 && s->right.lab[min + cf] == spec) {
			swap_labels(&s->right, min + cf, cb);
			min = orbit_prune(s);
		}
	}
	return min;
}

static void
rewind_coloring(struct saucy *s, struct coloring *c, int lev)
{
	int i, cf, ff, splits = s->splitlev[lev];
	for (i = s->nsplits - 1; i >= splits; --i) {
		cf = s->splitfrom[i];
		ff = s->splitwho[i];
		c->clen[cf] += c->clen[ff] + 1;
		fix_fronts(c, cf, ff);
	}
}

static int
do_backtrack(struct saucy *s)
{
	int i, cf, cb;

	/* Undo the splits up to this level */
	rewind_coloring(s, &s->right, s->lev);
	s->nsplits = s->splitlev[s->lev];

	/* Rewind diff information */
	for (i = s->ndiffs - 1; i >= s->difflev[s->lev]; --i) {
		s->diffmark[s->diffs[i]] = 0;
	}
	s->ndiffs = s->difflev[s->lev];
	s->nundiffs = s->undifflev[s->lev];

	/* Point to the target cell */
	cf = s->start[s->lev];
	cb = cf + s->right.clen[cf];

	/* Update ancestor with zeta if we've rewound more */
	if (s->anc > s->lev) {
		s->anc = s->lev;
		s->indmin = s->left.lab[cb];
		s->match = 1;
		note_anctar_reps(s);
	}

	/* Perform backtracking appropriate to our location */
	return s->lev == s->anc
		? backtrack_leftmost(s)
		: backtrack_other(s);
}

static int
backtrack_loop(struct saucy *s)
{
	int min;

	/* Backtrack as long as we're exhausting target cells */
	for (--s->lev; s->lev; --s->lev) {
		min = do_backtrack(s);
		if (min != -1) return min + s->start[s->lev];
	}
	return -1;
}

static int
backtrack(struct saucy *s)
{
	int min, old, tmp;
	old = s->nsplits;
	min = backtrack_loop(s);
	tmp = s->nsplits;
	s->nsplits = old;
	rewind_coloring(s, &s->left, s->lev+1);
	s->nsplits = tmp;

	return min;
}

static int
backtrack_bad(struct saucy *s)
{
	int min, old, tmp;
	old = s->lev;
	min = backtrack_loop(s);	
	tmp = s->nsplits;
	s->nsplits = s->splitlev[old];
	rewind_coloring(s, &s->left, s->lev+1);
	s->nsplits = tmp;

	return min;
}

static void
prepare_permutation(struct saucy *s)
{
	int i, k;
	for (i = 0; i < s->ndiffs; ++i) {
		k = s->right.unlab[s->diffs[i]];
		s->unsupp[i] = s->left.lab[k];
		s->gamma[s->left.lab[k]] = s->right.lab[k];
	}
}

static void
unprepare_permutation(struct saucy *s)
{
	int i;
	for (i = 0; i < s->ndiffs; ++i) {
		s->gamma[s->unsupp[i]] = s->unsupp[i];
	}
}

static int
do_search(struct saucy *s)
{
	int min;

	unprepare_permutation(s);

	/* Backtrack to the ancestor with zeta */
	if (s->lev > s->anc) s->lev = s->anc + 1;

	/* Perform additional backtracking */
	min = backtrack(s);

	/* Keep going while there are tree nodes to expand */
	while (s->lev) {

		/* Descend to a new leaf node */	
		if (descend(s, &s->right, s->start[s->lev], min)
				&& descend_left(s)) {

			/* Prepare permutation */
			prepare_permutation(s);

			/* If we found an automorphism, return it */
			if (s->is_automorphism(s)) {
				++s->stats->gens;
				s->stats->support += s->ndiffs;
				update_theta(s);
				return s->consumer(s->n, s->gamma,
					s->ndiffs, s->unsupp, s->arg);
			}
			else {
				unprepare_permutation(s);
			}
		}

		/* If we get here, something went wrong; backtrack */
		++s->stats->bads;
		min = backtrack_bad(s);
	}

	/* Normalize group size */
	while (s->stats->grpsize_base >= 10.0) {
		s->stats->grpsize_base /= 10;
		++s->stats->grpsize_exp;
	}
	return 0;
}

void
saucy_search(
	struct saucy *s,
	const struct saucy_graph *g,
	int directed,
	const int *colors,
	saucy_consumer *consumer,
	void *arg,
	struct saucy_stats *stats)
{
	int i, j, max = 0;

	/* Save client information */
	s->stats = stats;
	s->arg = arg;
	s->consumer = consumer;

	/* Save graph information */
	s->n = g->n;
	s->adj = g->adj;
	s->edg = g->edg;
	s->dadj = g->adj + g->n + 1;
	s->dedg = g->edg + g->e;

	/* Polymorphism */
	if (directed) {
		s->is_automorphism = is_directed_automorphism;
		s->ref_singleton = ref_singleton_directed;
		s->ref_nonsingle = ref_nonsingle_directed;
	}
	else {
		s->is_automorphism = is_undirected_automorphism;
		s->ref_singleton = ref_singleton_undirected;
		s->ref_nonsingle = ref_nonsingle_undirected;
	}

	/* Initialize scalars */
	s->indmin = 0;
	s->lev = s->anc = 1;
	s->ndiffs = s->nundiffs = s->ndiffnons = 0;

	/* The initial orbit partition is discrete */
	for (i = 0; i < s->n; ++i) {
		s->theta[i] = i;
	}

	/* The initial permutation is the identity */
	for (i = 0; i < s->n; ++i) {
		s->gamma[i] = i;
	}

	/* Initially every cell of theta has one element */
	for (i = 0; i < s->n; ++i) {
		s->thsize[i] = 1;
	}

	/* Every theta rep list is singleton */
	for (i = 0; i < s->n; ++i) {
		s->thprev[i] = s->thnext[i] = i;
	}

	/* We have no pairs yet */
	s->npairs = 0;
	for (i = 0; i < s->n; ++i) {
		s->unpairs[i] = -1;
	}

	/* Ensure no stray pointers in undiffnons, which is checked by removed_diffnon() */
	for (i = 0; i < s->n; ++i) {
		s->undiffnons[i] = -1;
	}

	/* Initialize stats */
	s->stats->grpsize_base = 1.0;
	s->stats->grpsize_exp = 0;
	s->stats->nodes = 1;
	s->stats->bads = s->stats->gens = s->stats->support = 0;

	/* Prepare for refinement */
	s->nninduce = s->nsinduce = 0;
	s->csize = 0;

	/* Count cell sizes */
	for (i = 0; i < s->n; ++i) {
		s->ccount[colors[i]]++;
		if (max < colors[i]) max = colors[i];
	}
	s->nsplits = max + 1;

	/* Build cell lengths */
	s->left.clen[0] = s->ccount[0] - 1;
	for (i = 0; i < max; ++i) {
		s->left.clen[s->ccount[i]] = s->ccount[i+1] - 1;
		s->ccount[i+1] += s->ccount[i];
	}

	/* Build the label array */
	for (i = 0; i < s->n; ++i) {
		set_label(&s->left, --s->ccount[colors[i]], i);
	}

	/* Clear out ccount */
	for (i = 0; i <= max; ++i) {
		s->ccount[i] = 0;
	}

	/* Update refinement stuff based on initial partition */
	for (i = 0; i < s->n; i += s->left.clen[i]+1) {
		add_induce(s, &s->left, i);
		fix_fronts(&s->left, i, i);
	}

	/* Prepare lists based on cell lengths */
	for (i = 0, j = -1; i < s->n; i += s->left.clen[i] + 1) {
		if (!s->left.clen[i]) continue;
		s->prevnon[i] = j;
		s->nextnon[j] = i;
		j = i;
	}

	/* Fix the end */
	s->prevnon[s->n] = j;
	s->nextnon[j] = s->n;

	/* Preprocessing after initial coloring */
	s->split = split_init;
	refine(s, &s->left);

	/* Descend along the leftmost branch and compute zeta */
	descend_leftmost(s);
	s->split = split_other;

	/* Our common ancestor with zeta is the current level */
	s->stats->levels = s->anc = s->lev;

	/* Copy over this data to our non-leftmost coloring */
	memcpy(s->right.lab, s->left.lab, s->n * sizeof(int));
	memcpy(s->right.unlab, s->left.unlab, s->n * sizeof(int));
	memcpy(s->right.clen, s->left.clen, s->n * sizeof(int));
	memcpy(s->right.cfront, s->left.cfront, s->n * sizeof(int));

	/* The reps are just the labels at this point */
	memcpy(s->threp, s->left.lab, s->n * sizeof(int));
	memcpy(s->thfront, s->left.unlab, s->n * sizeof(int));

	/* Keep running till we're out of automorphisms */
	while (do_search(s));
}

static int *ints(int n) { return malloc(n * sizeof(int)); }
static int *zeros(int n) { return calloc(n, sizeof(int)); }
static char *bits(int n) { return calloc(n, sizeof(char)); }

struct saucy *
saucy_alloc(int n)
{
	struct saucy *s = malloc(sizeof(struct saucy));
	if (s == NULL) return NULL;

	s->ninduce = ints(n);
	s->sinduce = ints(n);
	s->indmark = bits(n);
	s->left.cfront = zeros(n);
	s->left.clen = ints(n);
	s->right.cfront = zeros(n);
	s->right.clen = ints(n);
	s->stuff = bits(n+1);
	s->bucket = ints(n+2);
	s->count = ints(n+1);
	s->ccount = zeros(n);
	s->clist = ints(n);
	s->nextnon = ints(n+1) + 1;
	s->prevnon = ints(n+1);
	s->anctar = ints(n);
	s->start = ints(n);
	s->gamma = ints(n);
	s->junk = ints(n);
	s->theta = ints(n);
	s->thsize = ints(n);
	s->left.lab = ints(n);
	s->left.unlab = ints(n);
	s->right.lab = ints(n);
	s->right.unlab = ints(n);
	s->splitwho = ints(n);
	s->splitfrom = ints(n);
	s->splitlev = ints(n+1);
	s->unsupp = ints(n);
	s->conncnts = zeros(n);
	s->diffmark = bits(n);
	s->diffs = ints(n);
	s->difflev = ints(n);
	s->undifflev = ints(n);
	s->specmin = ints(n);
	s->thnext = ints(n);
	s->thprev = ints(n);
	s->threp = ints(n);
	s->thfront = ints(n);
	s->pairs = ints(n);
	s->unpairs = ints(n);
	s->diffnons = ints(n);
	s->undiffnons = ints(n);

	if (s->ninduce && s->sinduce && s->left.cfront && s->left.clen
		&& s->right.cfront && s->right.clen
		&& s->stuff && s->bucket && s->count && s->ccount
		&& s->clist && s->nextnon-1 && s->prevnon
		&& s->start && s->gamma && s->theta && s->left.unlab
		&& s->right.lab && s->right.unlab
		&& s->left.lab && s->splitwho && s->junk
		&& s->splitfrom && s->splitlev && s->thsize
		&& s->unsupp && s->conncnts && s->anctar
		&& s->diffmark && s->diffs && s->indmark
		&& s->thnext && s->thprev && s->threp && s->thfront
		&& s->pairs && s->unpairs && s->diffnons && s->undiffnons
		&& s->difflev && s->undifflev && s->specmin)
	{
		return s;
	}
	else {
		saucy_free(s);
		return NULL;
	}
}

void
saucy_free(struct saucy *s)
{
	free(s->undiffnons);
	free(s->diffnons);
	free(s->unpairs);
	free(s->pairs);
	free(s->thfront);
	free(s->threp);
	free(s->thnext);
	free(s->thprev);
	free(s->specmin);
	free(s->anctar);
	free(s->thsize);
	free(s->undifflev);
	free(s->difflev);
	free(s->diffs);
	free(s->diffmark);
	free(s->conncnts);
	free(s->unsupp);
	free(s->splitlev);
	free(s->splitfrom);
	free(s->splitwho);
	free(s->right.unlab);
	free(s->right.lab);
	free(s->left.unlab);
	free(s->left.lab);
	free(s->theta);
	free(s->junk);
	free(s->gamma);
	free(s->start);
	free(s->prevnon);
	free(s->nextnon-1);
	free(s->clist);
	free(s->ccount);
	free(s->count);
	free(s->bucket);
	free(s->stuff);
	free(s->right.clen);
	free(s->right.cfront);
	free(s->left.clen);
	free(s->left.cfront);
	free(s->indmark);
	free(s->sinduce);
	free(s->ninduce);
	free(s);
}
