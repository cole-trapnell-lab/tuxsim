#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#endif

#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <string>
#include <vector>

#include "common.h"
#include "scaffolds.h"
#include "gtf_tracking.h"
#include "hits.h"

using namespace std;

const char *short_options = "";

static struct option long_options[] = {
	{"inner-dist-mean",			required_argument,       0,          'm'},
	{"inner-dist-stddev",		required_argument,       0,          's'},
	{0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "tuxsim v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   tuxsim [options] <ref_transcripts.gtf> <output_sample.sam>\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-m/--frag-length-mean         the mean fragment length                         [ default:     200 ]\n");
	fprintf(stderr, "-s/--frag-length-std-dev      the standard deviation of the fragment length    [ default:     50 ]\n");
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
			case 'm':
				frag_length_mean = (uint32_t)parseInt(-200, "-m/--frag-length-mean arg must be at least -200", print_usage);
				break;
			case 's':
				frag_length_std_dev = (uint32_t)parseInt(0, "-s/--frag-length-std-dev arg must be at least 0", print_usage);
				break;
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	
	return 0;
}

struct ScaffoldSorter
{
	ScaffoldSorter(RefSequenceTable& _rt) : rt(_rt) {} 
	bool operator()(const Scaffold& lhs, const Scaffold& rhs)
	{
		const char* lhs_name = rt.get_name(lhs.ref_id());
		const char* rhs_name = rt.get_name(rhs.ref_id());
		int c = strcmp(lhs_name, rhs_name);
		if (c != 0)
		{
			return c < 0;
		}
		else
		{
			return lhs.left() < rhs.left();
		}
	}
	
	RefSequenceTable& rt;
};

void load_ref_rnas(FILE* ref_mRNA_file, 
				   RefSequenceTable& rt,
				   vector<Scaffold>& ref_mRNAs) 
{
	GList<GSeqData> ref_rnas;
	
	if (ref_mRNA_file)
	{
		//read_mRNAs(ref_mRNA_file, false, ref_rnas, ref_rnas, NULL, -1, false);
		read_mRNAs(ref_mRNA_file, ref_rnas);
	}
	
	// Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			char* name = GffObj::names->gseqs.getName(ref_rnas[j]->gseq_id);
			RefID ref_id = rt.get_id(name, NULL);
			for (int i = 0; i < ref_rnas[j]->mrnas_f.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_f[i]);
				vector<AugmentedCuffOp> ops;
				
				for (int e = 0; e < rna.exons.Count(); ++e)
				{
					GffExon& ex = *(rna.exons[e]);
					ops.push_back(AugmentedCuffOp(CUFF_MATCH, ex.start - 1, ex.end - ex.start + 1));
					if (e + 1 < rna.exons.Count())
					{
						GffExon& next_ex = *(rna.exons[e+1]);
						ops.push_back(AugmentedCuffOp(CUFF_INTRON, ex.end, next_ex.start - ex.end - 1));
					}
				}
				
				Scaffold ref_scaff(ref_id, CUFF_FWD, ops);
				if (rna.getID())
				{
					ref_scaff.annotated_trans_id(rna.getID());
				}
				if (rna.getGene())
					ref_scaff.annotated_gene_id(rna.getGene());
				char* short_name = rna.getAttr("gene_name");
				if (short_name)
				{
					ref_scaff.annotated_gene_name(short_name);
				}
				
				char* nearest_ref_match = rna.getAttr("nearest_ref");
				char* class_code = rna.getAttr("class_code");
				
				if (nearest_ref_match && class_code)
				{
					ref_scaff.nearest_ref_id(nearest_ref_match);
					ref_scaff.nearest_ref_classcode(*class_code);
				}
				
				char* protein_id = rna.getAttr("p_id");
				if (protein_id)
				{
					ref_scaff.annotated_protein_id(protein_id);
				}
				
				char* tss_id = rna.getAttr("tss_id");
				if (tss_id)
				{
					ref_scaff.annotated_tss_id(tss_id);
				}
				
				ref_mRNAs.push_back(ref_scaff); 
			}
			
			for (int i = 0; i < ref_rnas[j]->mrnas_r.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_r[i]);
				vector<AugmentedCuffOp> ops;
				
				for (int e = 0; e < rna.exons.Count(); ++e)
				{
					GffExon& ex = *(rna.exons[e]);
					ops.push_back(AugmentedCuffOp(CUFF_MATCH, ex.start - 1, ex.end - ex.start + 1));
					if (e + 1 < rna.exons.Count())
					{
						GffExon& next_ex = *(rna.exons[e+1]);
						ops.push_back(AugmentedCuffOp(CUFF_INTRON, ex.end, next_ex.start - ex.end - 1));
					}
				}
				
				Scaffold ref_scaff(ref_id, CUFF_REV, ops);
				if (rna.getID())
				{
					ref_scaff.annotated_trans_id(rna.getID());
				}
				if (rna.getGene())
					ref_scaff.annotated_gene_id(rna.getGene());
				char* short_name = rna.getAttr("gene_name");
				if (short_name)
				{
					ref_scaff.annotated_gene_name(short_name);
				}
				
				char* nearest_ref_match = rna.getAttr("nearest_ref");
				char* class_code = rna.getAttr("class_code");
				
				if (nearest_ref_match && class_code)
				{
					ref_scaff.nearest_ref_id(nearest_ref_match);
					ref_scaff.nearest_ref_classcode(*class_code);
				}
				
				char* protein_id = rna.getAttr("p_id");
				if (protein_id)
				{
					ref_scaff.annotated_protein_id(protein_id);
				}
				
				char* tss_id = rna.getAttr("tss_id");
				if (tss_id)
				{
					ref_scaff.annotated_tss_id(tss_id);
				}
				
				ref_mRNAs.push_back(ref_scaff); 
			}
		}
		ScaffoldSorter sorter(rt);
		sort(ref_mRNAs.begin(), ref_mRNAs.end(), sorter);
	}
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string ref_gtf_filename = argv[optind++];
    	
	// seed the random number generator 
	srand48(time(NULL));
	
	FILE* ref_gtf = NULL;
	if (ref_gtf_filename != "")
	{
		ref_gtf = fopen(ref_gtf_filename.c_str(), "r");
		if (!ref_gtf)
		{
			fprintf(stderr, "Error: cannot open GTF file %s for reading\n",
					ref_gtf_filename.c_str());
			exit(1);
		}
	}
	
	return 0;
}