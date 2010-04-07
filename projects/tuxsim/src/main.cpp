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

//#include <boost/math/distributions/normal.hpp> 
//using boost::math::normal;

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/random/normal_distribution.hpp>
//#include <boost/random/variate_generator.hpp>

#include "common.h"
#include "scaffolds.h"
#include "gtf_tracking.h"

#include "GStr.h"
#include "GFaSeqGet.h"

#include "hits.h"

using namespace std;
using namespace boost;

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
    fprintf(stderr, "Usage:   tuxsim [options] <ref_transcripts.gtf> <fasta_di> <out_sample_name>\n");
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

char* getGSeqName(int gseq_id)
{
	return GffObj::names->gseqs.getName(gseq_id);
}

char* getFastaFile(int gseq_id) 
{
	if (fastadir==NULL) return NULL;
	GStr s(fastadir);
	s.trimR('/');
	s.appendfmt("/%s",getGSeqName(gseq_id));
	GStr sbase=s;
	if (!fileExists(s.chars())) 
		s.append(".fa");
	if (!fileExists(s.chars())) s.append("sta");
	if (fileExists(s.chars())) return Gstrdup(s.chars());
	else
	{
		GMessage("Warning: cannot find genomic sequence file %s{.fa,.fasta}\n",sbase.chars());
		return NULL;
	}
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


//FIXME: needs refactoring
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
	
	int last_gseq_id = -1;
	GFaSeqGet* faseq = NULL;
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
				
				if (rna.gseq_id != last_gseq_id)
				{
					delete faseq;
					faseq = NULL;
					last_gseq_id = rna.gseq_id;
					char* sfile = getFastaFile(last_gseq_id);
					if (sfile != NULL) 
					{
						if (verbose)
							GMessage("Processing sequence from fasta file '%s'\n",sfile);
						faseq = new GFaSeqGet(sfile, false);
						faseq->loadall();
						GFREE(sfile);
					}
					else 
					{
						assert (false);
					}
				}
				
				assert (faseq);
				
				int seqlen = 0;
				char* rna_seq = rna.getSpliced(faseq, false, &seqlen);
				
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
				
				ref_scaff.seq(rna_seq);
				GFREE(rna_seq);
				
				ref_mRNAs.push_back(ref_scaff); 
				
			}
			
			for (int i = 0; i < ref_rnas[j]->mrnas_r.Count(); ++i)
			{	
				GffObj& rna = *(ref_rnas[j]->mrnas_r[i]);
				
				if (rna.gseq_id != last_gseq_id)
				{
					delete faseq;
					faseq = NULL;
					last_gseq_id = rna.gseq_id;
					char* sfile = getFastaFile(last_gseq_id);
					if (sfile != NULL) 
					{
						if (verbose)
							GMessage("Processing sequence from fasta file '%s'\n",sfile);
						faseq = new GFaSeqGet(sfile, false);
						faseq->loadall();
						GFREE(sfile);
					}
					else 
					{
						assert (false);
					}
				}
				
				assert (faseq);
				
				int seqlen = 0;
				char* rna_seq = rna.getSpliced(faseq, false, &seqlen);
				
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
				
				ref_scaff.seq(rna_seq);
				GFREE(rna_seq);
				
				ref_mRNAs.push_back(ref_scaff); 
			}
		}
		ScaffoldSorter sorter(rt);
		sort(ref_mRNAs.begin(), ref_mRNAs.end(), sorter);
	}
}

void assign_expression_ranks(const vector<Scaffold>& ref_mRNAs,
							 vector<unsigned int>& expr_rank)
{
	for (size_t i = 0; i < expr_rank.size(); ++i)
	{
		expr_rank[i] = (unsigned int)i;
	}
	
	random_shuffle(expr_rank.begin(), expr_rank.end());
}


struct FluxRankAbundancePolicy
{
	FluxRankAbundancePolicy(double x_0, double k, double x_1) : 
		_x_0(x_0), 
		_k(k), 
		_x_1(x_1) {}
	
	double rho(unsigned int rank) const
	{
		double x = rank + 1;
		double e = exp((x / _x_1) * (-1 - (x / _x_1))); 
		double r = pow((x / _x_0), _k) * e;
		return r;
	} 
	
private:
	double _x_0;
	double _k;
	double _x_1;
};

// Assigns the transcriptome wide relative abundances
template<class RankedBasedAbundancePolicy>
void assign_abundances(const vector<Scaffold>& ref_mRNAs,
					   const RankedBasedAbundancePolicy& rank_based,
					   vector<double>& expr_rho)
{
	vector<unsigned int> expr_rank(ref_mRNAs.size(), 0);
	
	assert(expr_rho.size() == expr_rank.size());
	
	assign_expression_ranks(ref_mRNAs, expr_rank);
	
	double total_mols = 0.0;
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		int rank = expr_rank[i];
		double rho = rank_based.rho(rank);
		expr_rho[i] = rho;
		total_mols += rho;
	}
	
	assert (total_mols > 0.0);
	
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		expr_rho[i] /= total_mols;
	}
}

// FIXME: should use effective length
void calc_frag_abundances(const vector<Scaffold>& ref_mRNAs,
						  const vector<double>& expr_rho,
						  vector<double>& expr_alpha)
{
	double fragment_pool_size = 0.0;
	
	assert (expr_alpha.size() == ref_mRNAs.size() &&
			expr_alpha.size() == expr_rho.size());
	
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		fragment_pool_size += ref_mRNAs[i].length() * expr_rho[i];
	}
	
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		double frag_n = ref_mRNAs[i].length() * expr_rho[i];
		expr_alpha[i] = frag_n / fragment_pool_size;
	}
}

struct FastqOutfilePair
{
	FastqOutfilePair(FILE* left, FILE* right) : 
	left_out_file(left), 
	right_out_file(right) {}
	
	FILE* left_out_file;
	FILE* right_out_file;
};

struct LibraryFragment
{
	const Scaffold* source_seq;
	
	// Fragment specified as an open interval [start, end)
	int start;
	int end;	
};

class FragmentPolicy
{
public:
	virtual void next_fragment(const Scaffold& molecule,
							   LibraryFragment& fragment) = 0;
};

class NormalFragments : public FragmentPolicy
{
	
	// This is a typedef for a random number generator.
	// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
	typedef boost::minstd_rand base_generator_type;
	
	typedef variate_generator<base_generator_type&, boost::uniform_real<> > uniform_generator_type;
	typedef variate_generator<base_generator_type&, normal_distribution<> > normal_generator_type;
	
	
	//boost::uniform_real<> uni_dist(0,1);
	// uni(generator, uni_dist);
	
public:
	
	NormalFragments(int mean_frag_length, int frag_length_sd) :
		_base_generator(base_generator_type(time(NULL))),
		_uniform_generator(uniform_generator_type(_base_generator, uniform_real<>(0,1))),
		_length_generator(normal_generator_type(_base_generator, 
												normal_distribution<>(mean_frag_length, frag_length_sd)))
	{
	}
	
	virtual void next_fragment(const Scaffold& molecule,
							   LibraryFragment& fragment)
	{
		int frag_length = _length_generator();
		int frag_start = (molecule.length() - frag_length) * 
			_uniform_generator();
		assert (frag_start < molecule.length());
		
		int frag_end = frag_start + frag_length;
		assert (frag_end < molecule.length());
		
		fragment.source_seq = &molecule;
		fragment.start = frag_start;
		fragment.end = frag_end;
	}
	
private:
	base_generator_type _base_generator;
	uniform_generator_type _uniform_generator;
	normal_generator_type _length_generator;
};

typedef vector<shared_ptr<ReadHit> > ReadsForFragment; 

class SequencingPolicy
{
public:
	virtual bool reads_for_fragment(const LibraryFragment& frag, 
									ReadsForFragment& reads) = 0;
};

/*******************************************************************************
                            Sequencing Policies
*******************************************************************************/
struct IlluminaChIPSeqPE : public SequencingPolicy
{
	
	IlluminaChIPSeqPE(int left_read_len, 
					  int right_read_len) :
		_left_len(left_read_len),
		_right_len(right_read_len) {}
	
	bool reads_for_fragment(const LibraryFragment& frag, 
							ReadsForFragment& reads)
	{
		int frag_length = frag.end - frag.start;
		if (frag_length < _left_len ||
			frag_length < _right_len)
		{
			return false;
		}
		
		// Generate a new read
		
		return true;
	}
	
private:
	int _left_len;
	int _right_len;
};
/******************************************************************************/

class AssayProtocol
{
public:
	AssayProtocol(FragmentPolicy& frag_impl, 
				  SequencingPolicy& seq_impl) :
		_frag_impl(frag_impl), 
		_seq_impl(seq_impl) {}
	
	void next_reads(const Scaffold& molecule, ReadsForFragment& reads)
	{
		LibraryFragment frag;
		_frag_impl.next_fragment(molecule, frag);
		_seq_impl.reads_for_fragment(frag, reads);
	}
	
private:
	FragmentPolicy& _frag_impl;
	SequencingPolicy& _seq_impl;
};

void generate_reads(const vector<Scaffold>& ref_mRNAs, 
					const vector<double>& frag_abundances,
					int total_frags,
					FILE* sam_frag_out,
					FastqOutfilePair& fastq_out)
{
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
		int num_frags_for_mRNA = frag_abundances[i] * total_frags;
		for (size_t j = 0; j < num_frags_for_mRNA; ++j)
		{
			
		}
	}
}

void driver(FILE* ref_gtf, 
			FILE* sam_out,
			FastqOutfilePair& fastq_out)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	vector<Scaffold> ref_mRNAs;
	load_ref_rnas(ref_gtf, rt, ref_mRNAs);
	if (ref_mRNAs.empty())
		return;
	
	FluxRankAbundancePolicy flux_policy(5e7, -0.6, 9500);
	
	vector<double> expr_rho(ref_mRNAs.size(), 0.0);
	
	assign_abundances(ref_mRNAs,
					  flux_policy,
					  expr_rho);
	
	vector<double> expr_alpha(ref_mRNAs.size(), 0.0);
	calc_frag_abundances(ref_mRNAs,
						 expr_rho,
						 expr_alpha);
	
	int total_frags = 100;
	
	generate_reads(ref_mRNAs,
				   expr_alpha,
				   total_frags,
				   sam_out,
				   fastq_out);
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

	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    fastadir = argv[optind++];
	
	if(optind >= argc)
    {
        print_usage();
        return 1;
    }
	
    string out_prefix = argv[optind++];
    	
	// seed the random number generator 
	srand(time(NULL));
	
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
	
	FILE* sam_out = NULL;
	string out_sam_filename = out_prefix + ".sam";
	sam_out = fopen(out_sam_filename.c_str(), "w");
	if (!sam_out)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
				out_sam_filename.c_str());
		exit(1);
	}

	FILE* left_fastq_out = NULL;
	string left_out_fastq_filename = out_prefix + "_1.fq";
	left_fastq_out = fopen(left_out_fastq_filename.c_str(), "w");
	if (!sam_out)
	{
		fprintf(stderr, "Error: cannot open FASTQ file %s for writing\n",
				left_out_fastq_filename.c_str());
		exit(1);
	}
	
	FILE* right_fastq_out = NULL;
	string right_out_fastq_filename = out_prefix + "_1.fq";
	right_fastq_out = fopen(right_out_fastq_filename.c_str(), "w");
	if (!sam_out)
	{
		fprintf(stderr, "Error: cannot open FASTQ file %s for writing\n",
				right_out_fastq_filename.c_str());
		exit(1);
	}
	
	FastqOutfilePair fastq_out(left_fastq_out, right_fastq_out);
	
	driver(ref_gtf, sam_out, fastq_out);
	
	return 0;
}