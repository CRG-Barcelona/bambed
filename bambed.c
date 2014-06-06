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
  "   -shift=n                      shift read n bases\n"
  "   -length=l                     override insert-size length calculation from bam with l\n"
  "   -mapq=q                       minimum mapping-quality read to consider (default 20)\n"
  "   -include-duplicates=d         use up to d duplicated reads/fragments/pairs\n"
  "                                 using this option relies on bams being tagged with the\n"
  "                                 ZD tag (default uses bam flag for duplicates)\n"
  "   -include-B-reads              include \"B-reads\", tagged \"ZL\" in the bam\n"
  "   -include-bad-regions          include reads tagged \"ZR\" in the bam\n"
  "   -rg-whitelist=rg1,rg2,...     perform calculation using specific read-groups in the bam.\n"
  "   -rg-blacklist=rg1,rg2,...     exclude specfic read-groups (not compatible with -rg-whitelist)\n"
  "   -flag-counts=file             (for debugging) output counts of flags of reads used to file\n"              
  "   -count                        output a bed4 with counts for each region specified\n"
  "   -cc                         output split read bed for chrom capture bams (3C, Hi-C, etc)\n"
  "   -cc-inter                   output only inter-chromosomal pairs\n"
  "   -cc-intra                   output only intra-chromosomal pairs\n"
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
   {"cc", OPTION_BOOLEAN},
   {"cc-inter", OPTION_BOOLEAN},
   {"cc-intra", OPTION_BOOLEAN},
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

enum hicout
{
    both = 0,
    inter = 1,
    intra = 2,
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

void output_hic(struct metaBig *mb, char *outputfile)
/* updated to use pairbed */
{
    FILE *output = mustOpen(outputfile, "w");
    enum hicout h = both;
    struct bed *section;
    if (optionExists("cc-inter"))
	h = inter;
    else if (optionExists("cc-intra"))
	h = intra;    
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct lm *lm = lmInit(0);
	struct pairbed *pbList = metaBigPairbedFetch(mb, section->chrom, section->chromStart, section->chromEnd, lm);
	struct pairbed *pb;
	for (pb = pbList; pb != NULL; pb = pb->next)
	{
	    if ((h == both) || ((h == inter) && (!sameString(pb->chrom, pb->mChrom))) || 
		((h == intra) && (sameString(pb->chrom, pb->mChrom))))
	    {
		fprintf(output, "%s\t%d\t%s\t%s\t%d\t%s\n", pb->chrom, pb->chromStart, pb->strand, 
			pb->mChrom, pb->mChromStart, pb->mstrand);
	    }
	}
	lmCleanup(&lm);
    }
    carefulClose(&output);
}

void bambed(char *bigfile, char *outputfile)
/* bambed - main */
{
    char *regionsBed = optionVal("regions", NULL);
    boolean verbose = optionExists("verbose");
    if (verbose)
	printf("loading %s... \n", bigfile);
    struct metaBig *mb = metaBigOpen(bigfile, regionsBed);
    if (!mb)
	errAbort("could not load %s for some reason", bigfile);
    if (verbose)
	printf("loaded %s ok\n", bigfile);
    boolean do_cc = FALSE;
    if (optionExists("cc") || optionExists("cc-inter") || optionExists("cc-intra"))
	do_cc = TRUE;
    if (optionExists("cc-inter") && optionExists("cc-intra"))
	errAbort("-cc-inter and -cc-intra can't both be used");
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
    if (do_cc)
	output_hic(mb, outputfile);
    else if (optionExists("count"))
	output_counts(mb, outputfile);
    else
	output_bed(mb, outputfile);
    metaBigClose(&mb);
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
