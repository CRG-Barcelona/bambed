/* bambed - make a vector for use with those typical +/- plots */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include <beato/metaBig.h>
#include <sam.h>
#include <bam.h>

void usage()
/* Explain usage and exit. */
{
errAbort(
  "bambed - Convert BAM regions to BED6\n"
  "usage:\n"
  "   bambed input.bam:chr:start-end output.bed\n"
  "options:\n"
  "   -name=[dup|seq|fqh|qual]      name the bed items by sequence or original\n"
  "                                 header from the fastq, or quality\n"
  "   -regions=bed                  provide a bed of specific regions to extract from\n"
  "                                 the bam\n"
  "   -regions-every=size           Instead of specifying a bed of regions, specify a size\n"
  "                                 for regularly-spaced intervals.\n"
/*   "   -fifty                        try to improve the correct retrieval by taking only ones\n" */
/*   "                                 that are > 50%% with a region covering the start, and >= 50%%\n" */
/*   "                                 within a region on the end.\n" */
/* default */
  "   -shift=n                      shift read n bases\n"
  "   -length=l                     override insert-size length calculation from bam with l\n"
  "   -mapq=q                       minimum mapping-quality read to consider (default 20)\n"
  "   -include-duplicates=d         use up to d duplicated reads/fragments/pairs\n"
  "                                 using this option relies on bams being tagged with the\n"
  "                                 ZD tag (default uses bam flag for duplicates)\n"
  "   -include-B-reads              include \"B-reads\", tagged \"ZL\" in the bam\n"
  "   -include-bad-regions          include reads tagged \"ZR\" in the bam\n"
/*   "   -lift-artifical               use regions from -regions to change coordinates to\n" */
/*   "                                 be an artificial chromosome for testing\n" */
/*   "   -lift-gap=n                   if using -lift-artificial, make the buffer between\n" */
/*   "                                 regions n bases long (default 50)\n" */
  "   -rg-whitelist=rg1,rg2,...     perform calculation using specific read-groups in the bam.\n"
  "   -rg-blacklist=rg1,rg2,...     exclude specfic read-groups (not compatible with -rg-whitelist)\n"
  "   -flag-counts=file             (for debugging) output counts of flags of reads used to file\n"              
  "   -count                        output a bed4 with counts for each region specified\n"
  "   -hi-c                         output split read bed for Hi-C\n"
/*   "   -pe-unpaired                  unpair the paired-end reads (requires -unextended)\n" */
  "   -stranded                     for paired-end reads that are stranded i.e. the first one in\n"
  "                                 the pair should have the strand.\n"
  "   -verbose                      print some progress info\n"
  );
}

static struct optionSpec options[] = 
{
   {"name-type", OPTION_STRING},
   {"shift", OPTION_INT},
   {"length", OPTION_INT},
   {"mapq", OPTION_INT},
   {"lift-artificial", OPTION_BOOLEAN},
   {"lift-gap", OPTION_INT},
   {"strand", OPTION_STRING},
   {"count", OPTION_BOOLEAN},
   {"include-duplicates", OPTION_INT},
   {"include-B-reads", OPTION_BOOLEAN},
   {"include-bad-regions", OPTION_BOOLEAN},
   {"flag-counts", OPTION_STRING},
   {"hi-c", OPTION_BOOLEAN},
   {"regions", OPTION_STRING},
   {"regions-every", OPTION_INT},
   {"rg-whitelist", OPTION_STRING},
   {"rg-blacklist", OPTION_STRING},
   {"favorites", OPTION_STRING},
   {"stranded", OPTION_BOOLEAN},
   {"verbose", OPTION_BOOLEAN},
   {"fifty", OPTION_BOOLEAN},
   {NULL, 0},
};

void output_bed(struct metaBig *mb, char *outputfile)
/* just output it */
{
    FILE *output = mustOpen(outputfile, "w");
    boolean do_sort = optionExists("length") || optionExists("shift");
    struct bed *section;
    char *flag_counts_file = optionVal("flag-counts", NULL);
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct lm *lm = lmInit(0);
	struct bed6 *beds = metaBigBed6Fetch(mb, section->chrom, section->chromStart, section->chromEnd, lm);
	struct bed6 *bed;
	int size = section->chromEnd - section->chromStart;
	if (do_sort)
	    slSort(&beds, bed6Cmp);
	for (bed = beds; bed != NULL; bed = bed->next)
	    fprintf(output, "%s\t%d\t%d\t%s\t%d\t%c\n", bed->chrom, bed->chromStart, bed->chromEnd, bed->name, bed->score, bed->strand[0]);
	lmCleanup(&lm);
    }
    if (flag_counts_file)
	metaBigPrintFlagCounts(mb, flag_counts_file, FALSE);
    carefulClose(&output);
}

void output_counts(struct metaBig *mb, char *outputfile)
/* only counting */
{
    boolean fifty = optionExists("fifty");
    FILE *output = mustOpen(outputfile, "w");
    struct bed *section;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	long count = metaBigFiftyCount(mb, section->chrom, section->chromStart, section->chromEnd);
	fprintf(output, "%s\t%d\t%d\t%ld\n", section->chrom, section->chromStart, section->chromEnd, count);
    }
    carefulClose(&output);
}

void hiCBam(char *bigfile, char *outputfile)
{
    FILE *out = mustOpen(outputfile, "w");
    bam1_t *b[2];
    bam1_t *cur = NULL;
    bam1_t *pre = NULL;
    int curr = 0; 
    boolean has_prev = FALSE;
    if (!bigfile)
	errAbort("no file");
    bamFile in = (strcmp(bigfile, "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(bigfile, "r");
    bam_header_t *header = bam_header_read(in);
    b[0] = bam_init1();
    b[1] = bam_init1();
    while (bam_read1(in, b[curr]) >= 0)
    {
	cur = b[curr];
	pre = b[1-curr];
	if (has_prev)
	{
	    char *name1 = bam1_qname(cur);
	    char *name2 = bam1_qname(pre);
	    if (name1 && name2 && sameString(name1, name2))
	    /* we have both in a pair of paired-end reads */
	    {
		bam1_core_t *core_f = &pre->core;
		bam1_core_t *core_s = &cur->core;
		char *chrom_f = header->target_name[core_f->tid];
		char *chrom_s = header->target_name[core_s->tid];
		int start_f = core_f->pos;
		int start_s = core_s->pos;
		int end_f = bam_calend(core_f, bam1_cigar(pre));
		int end_s = bam_calend(core_s, bam1_cigar(cur));
		char strand_f = '+';
		char strand_s = '+';
		if (core_f->flag & BAM_FREVERSE)
		    strand_f = '-';
		if (core_s->flag & BAM_FREVERSE)
		    strand_s = '-';
		fprintf(out, "%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n", chrom_f, start_f, end_f, strand_f, chrom_s, start_s, end_s, strand_s);
		has_prev = FALSE;
	    }
	}
	else 
	    has_prev = TRUE;
	curr = 1 - curr;
    }
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    bam_close(in);
    carefulClose(&out);
}

void bambed(char *bigfile, char *outputfile)
/* bambed - main */
{
    char *regionsBed = optionVal("regions", NULL);
    boolean verbose = optionExists("verbose");
    if (bigfile == NULL)
	uglyf("WTF\n");
    if (optionExists("hi-c"))
	hiCBam(bigfile, outputfile);
    else
    {
	if (verbose)
	    printf("loading %s... \n", bigfile);
	struct metaBig *mb = metaBigOpen(bigfile, regionsBed);
	if (verbose)
	    printf("loaded %s ok\n", bigfile);
	if (optionExists("regions-every"))
	{
	    if (regionsBed)
		errAbort("cannot specify -regions-every with -regions");
	    else
		mb->sections = metaBig_chopGenome(mb, optionInt("regions-every", 10000));
	}	    
	enum metaBigNameType name_type = justADot;
	char *blacklist = optionVal("rg-blacklist", NULL);
	char *whitelist = optionVal("rg-whitelist", NULL);
	if (whitelist && blacklist)
	    errAbort("cannot use -rg-whitelist with -rg-blacklist");
	if (blacklist)
	    metaBigSetRgList(mb, blacklist, TRUE);
	else if (whitelist)
	    metaBigSetRgList(mb, whitelist, FALSE);
	int shift = optionInt("shift", 0);
	int length = optionInt("length", 0);
	char strand = 0;
	char *name_type_s = optionVal("name-type", NULL);
	if (name_type_s)
	{
	    if (sameWord(name_type_s, "seq"))
		name_type = sequence;
	    else if (sameWord(name_type_s, "fqh"))
		name_type = basicName;
	    else if (sameWord(name_type_s, "qual"))
		name_type = quality;
	    else if (sameWord(name_type_s, "dup"))
		name_type = duplicates;
	}
	if (optionExists("strand"))
	{
	    char *s = optionVal("strand", NULL);
	    if (s && (s[0] == '+'))
		strand = '+';
	    else
		strand = '-';
	}
	metaBigSetPositionalOptions(mb, length, shift, strand);
	metaBigSetNameOption(mb, name_type);
	mb->includeB = optionExists("include-B-reads");
	mb->includeBadRegions = optionExists("include-bad-regions");
	mb->useDupes = optionInt("include-duplicates", 0);
	mb->mapQ = optionInt("mapq", 20);
	mb->strandedPe = optionExists("stranded");
	if (optionExists("mapq"))
	    mb->useMapQ = TRUE;
	if (optionExists("count"))
	    output_counts(mb, outputfile);
	else
	    output_bed(mb, outputfile);
	metaBigClose(&mb);
    }
}

int main(int argc, char *argv[])
/* Process command line. */
{
optionInit(&argc, argv, options);
if (argc != 3)
    usage();
bambed(argv[1], argv[2]);
return 0;
}
