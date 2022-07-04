#ifndef SAUCY_H
#define SAUCY_H

#define SAUCY_VERSION "2.0"

typedef int saucy_consumer(int, const int *, int, int *, void *);

struct saucy;

struct saucy_stats {
	double grpsize_base;
	int grpsize_exp;
	int levels;
	int nodes;
	int bads;
	int gens;
	int support;
};

struct saucy_graph {
	int n;
	int e;
	int *adj;
	int *edg;
};

struct saucy *saucy_alloc(int n);

void saucy_search(
	struct saucy *s,
	const struct saucy_graph *graph,
	int directed,
	const int *colors,
	saucy_consumer *consumer,
	void *arg,
	struct saucy_stats *stats);

void saucy_free(struct saucy *s);

#endif
