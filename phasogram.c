#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <beato/bigs.h>
#include <beato/metaBig.h>

#include "bambed.h"

struct phase
/* all we need to keep track of for the phasogram construction */
{
    struct phase *next;
    int size;
    int *dist_hist;
};

static struct phase *alloc_phase(int size)
/* pog constructor, essentially */
{
    struct phase *pog;
    if (size < 2)
	errAbort("phasogram is too small.  make -max-phase bigger");
    AllocVar(pog);
    pog->size = size + 1;
    AllocArray(pog->dist_hist, pog->size);
    return pog;
}

static void phaseFree(struct phase **pPog)
/* not so exciting */
{
    struct phase *pog = *pPog;
    if (pog->dist_hist)
	freeMem(pog->dist_hist);
    freez(pPog);
}

static void output_phasogram(FILE *f, struct phase *pog)
/* output tab-formatted phasogram */
{
    int i;
    for (i = 1; i < pog->size; i++)
	fprintf(f, "%d\t%d\n", i, pog->dist_hist[i]);
}

/*** the real program ****/

static void add_to_hists(struct phase *pog, int *starts, int *counts, int array_len, int max_phase)
/* important function.  does most of the algorithm.  does the triangular looping */
/* through the array, making the counts.  */
{
    int i, j;
    for (i = 0; i < array_len - 1; i++)
    {
	if ((array_len < 2) || (starts[array_len-1] - starts[i] < max_phase))
	    break;
	for (j = i + 1; j < array_len; j++)
	{
	    int d = starts[j] - starts[i];
	    if (counts[i] < 1)
		break;
	    if ((d <= max_phase) && (counts[j] >= 1))
		pog->dist_hist[d] += 1.0;
	    else if (d > max_phase)
		break;
	}
    }
}

static void add_to_phasogram(struct phase *total, struct phase *addition)
/* add phasogram to another. for example to combine all chroms */
{
    int i;
    for (i = 1; i < total->size; i++)
	total->dist_hist[i] += addition->dist_hist[i];
}

static struct phase *build_phasogram(struct starts *starts, int max_phase)
/* the general algorithm of building a phasogram. */
{
    struct phase *pog = alloc_phase(max_phase);
    /* first make counts for histogram */
    add_to_hists(pog, starts->pos_starts, starts->pos_counts, starts->num_pos_starts, max_phase);
    add_to_hists(pog, starts->neg_starts, starts->neg_counts, starts->num_neg_starts, max_phase);
    return pog;
}

struct phase *build_pe_phasogram(struct metaBig *mb, int max_phase)
{
    struct phase *total_phase = alloc_phase(max_phase);
    struct bed *bed;
    int subset = optionInt("stack-subset", 0);
    for (bed = mb->sections; bed != NULL; bed = bed->next)
    {
	struct middles *mids = metaBig_get_middles(mb, bed->chrom, bed->chromStart, bed->chromEnd);
	if (mids)
	{
	    struct phase *pog = alloc_phase(max_phase);
	    add_to_hists(pog, mids->mids, mids->counts, mids->num_mids, max_phase);
	    add_to_phasogram(total_phase, pog);
	    phaseFree(&pog);
	    middlesFreeList(&mids);
	}
    }
    return total_phase;
}

static struct phase *build_se_phasogram(struct metaBig *mb, int max_phase)
{
    struct phase *total_phase = alloc_phase(max_phase);
    struct bed *bed;
    for (bed = mb->sections; bed != NULL; bed = bed->next)
    {
	struct starts *starts = metaBig_get_starts(mb, bed->chrom, bed->chromStart, bed->chromEnd);
	if (starts)
	{
	    struct phase *pog = build_phasogram(starts, max_phase);
	    add_to_phasogram(total_phase, pog);
	    phaseFree(&pog);
	    startsFreeList(&starts);
	}
    }
    return total_phase;
}

void do_phasogram(struct metaBig *mb, char *outputfile, int max_phase)
/* the main part of the program.  loop through the input, */
/* do something, write the output */
{
    struct phase *total_phase = alloc_phase(max_phase);
    FILE *out = mustOpen(outputfile, "w");
    if (mb->pe)
	total_phase = build_pe_phasogram(mb, max_phase);
    else
	total_phase = build_se_phasogram(mb, max_phase);
    output_phasogram(out, total_phase);
    phaseFree(&total_phase);
    carefulClose(&out);
}
