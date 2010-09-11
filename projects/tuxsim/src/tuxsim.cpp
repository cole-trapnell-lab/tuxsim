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

#include <iostream>
#include <locale>


//#include <boost/math/distributions/normal.hpp> 
//using boost::math::normal;

#include "common.h"
#include "scaffolds.h"
#include "gtf_tracking.h"

#include "GStr.h"
#include "GFaSeqGet.h"
#include "GFastaFile.h"

#include "options.h"

#include "hits.h"
#include "abundances.h"
#include "fragments.h"
#include "sequencing.h"

using namespace std;
using namespace boost;

char* getGSeqName(int gseq_id)
{
	return GffObj::names->gseqs.getName(gseq_id);
}

char* getFastaFile(int gseq_id) 
{
	if (fastadir == "") return NULL;
	GStr s(fastadir.c_str());
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
        uint32_t l_id = lhs.ref_id();
        uint32_t r_id = rhs.ref_id();
        
        //uint32_t l_len = rt.get_len(lhs.ref_id());
        //uint32_t r_len = rt.get_len(rhs.ref_id());
        if (l_id != r_id)
        {
            //if (l_len != 0 && r_len != 0)
            //    return l_len > r_len;
            //else 
            return (strcmp(rt.get_name(lhs.ref_id()), rt.get_name(rhs.ref_id())) < 0);
        }
        return lhs.left() < rhs.left();
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
    
    sort(ref_mRNAs.begin(), ref_mRNAs.end(), ScaffoldSorter(rt));
	
	int last_gseq_id = -1;
	GFaSeqGet* faseq = NULL;
	// Geo groups them by chr.
	if (ref_rnas.Count()>0) //if any ref data was loaded
	{
		for (int j = 0; j < ref_rnas.Count(); ++j) 
		{    //ref data is grouped by genomic sequence
			char* name = GffObj::names->gseqs.getName(ref_rnas[j]->gseq_id);
			
			RefID ref_id = rt.get_id(name, NULL, 0);
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
                string rs = rna_seq;
				std::transform(rs.begin(), rs.end(), rs.begin(), (int (*)(int))std::toupper);
				ref_scaff.seq(rs);
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
				string rs = rna_seq;
                std::transform(rs.begin(), rs.end(), rs.begin(), (int (*)(int))std::toupper);
                reverse_complement(rs);
				ref_scaff.seq(rs);
				GFREE(rna_seq);
				
				ref_mRNAs.push_back(ref_scaff); 
			}
		}
		ScaffoldSorter sorter(rt);
		sort(ref_mRNAs.begin(), ref_mRNAs.end(), sorter);
	}
	delete faseq;
}

struct FastqOutfilePair
{
	FastqOutfilePair(FILE* left, FILE* right) : 
	left_out_file(left), 
	right_out_file(right) {}
	
	FILE* left_out_file;
	FILE* right_out_file;
};


/******************************************************************************/

class AssayProtocol
{
public:
	AssayProtocol(FragmentPolicy& frag_impl, 
				  SequencingPolicy& seq_impl) :
		_frag_impl(frag_impl), 
		_seq_impl(seq_impl) {}
	
	bool next_reads(const Scaffold& molecule, ReadsForFragment& reads)
	{
		LibraryFragment frag;
		
		if (_frag_impl.next_fragment(molecule, frag))
		{
			return _seq_impl.reads_for_fragment(frag, reads);
		}
		return false;
	}
    
    const FragmentPolicy& fragment_policy() const { return _frag_impl; }
    const SequencingPolicy& sequencing_policy() const { return _seq_impl; }
	
private:
	FragmentPolicy& _frag_impl;
	SequencingPolicy& _seq_impl;
};

struct SortReads
{
    bool operator()(shared_ptr<ReadHit> lhs, shared_ptr<ReadHit> rhs)
    {
        return lhs->left() < rhs->left();
    }
};

void print_sam_header(FILE* sam_out, const RefSequenceTable& rt)
{
    fprintf(sam_out, "@HD\tVN:1.0\tSO:sorted\n");
    map<string, uint32_t> seq_dict;
    for (RefSequenceTable::const_iterator itr = rt.begin();
         itr != rt.end();
         ++itr)
    {
//        fprintf(sam_out, 
//                "@SQ\tSN:%s\tLN:%u\n", 
//                itr->second.name, 
//                itr->second.len);
        seq_dict[itr->second.name] = itr->second.len;
    }
    
    for (map<string, uint32_t>::const_iterator itr = seq_dict.begin();
         itr != seq_dict.end();
         ++itr)
    {
        fprintf(sam_out, 
                "@SQ\tSN:%s\tLN:%u\n", 
                itr->first.c_str(), 
                itr->second);
    }
    
    fprintf(sam_out, "@PG\tID:TuxSim\tVN:%s", PACKAGE_VERSION);
}

void print_fastq_read(const ReadHit& read,
                      FILE* fastq_file)
{
    fprintf(fastq_file, "@%s\n",read.name().c_str());
    fprintf(fastq_file, "%s\n", read.seq().c_str());
    fprintf(fastq_file, "+\n");
    fprintf(fastq_file, "%s\n", read.qual().c_str());
}

void print_fastq_pair(const ReadsForFragment& reads, 
                      FastqOutfilePair& fastq_out)
{
    assert (reads.size() == 2);
	
	if (reads.front()->sam_flag() & BAM_FREAD1)
	{
		print_fastq_read(*(reads.front()), fastq_out.left_out_file);
		print_fastq_read(*(reads.back()), fastq_out.right_out_file);
	}
	else
	{
		print_fastq_read(*(reads.back()), fastq_out.left_out_file);
		print_fastq_read(*(reads.front()), fastq_out.right_out_file);
    }
}

void print_aligned_read(const ReadHit& read,
                        RefSequenceTable& rt,
                        FILE* sam_out)
{
    const char* read_name = read.name().c_str();
    
    int sam_flag = read.sam_flag();
    const char* ref_name = rt.get_name(read.ref_id());
    int ref_pos = read.left();
    const vector<CigarOp>& cigar = read.cigar();
    string cigar_str;
    for (size_t i = 0; i < cigar.size(); ++i)
    {
        char buf[256];
        char opcode = 0;
        switch(cigar[i].opcode)
        {
            case MATCH:
                opcode = 'M';
                break;
            case INS:
                opcode = 'I';
                break;
            case DEL:
                opcode = 'D';
                break;
            case REF_SKIP:
                opcode = 'N';
                break;
            case SOFT_CLIP:
                opcode = 'S';
                break;
            case HARD_CLIP:
                opcode = 'H';
                break;
            case PAD:
                opcode = 'P';
                break;
        }
        sprintf(buf, "%d%c", cigar[i].length, opcode);
        cigar_str += buf;
    }
    
    const char* mate_ref_name = rt.get_name(read.partner_ref_id());
    int mate_ref_pos = read.partner_pos();
    
    string seq = read.seq();
    string qual = read.qual();
    
    if (sam_flag & BAM_FREVERSE)
    {
        reverse_complement(seq);
        reverse(qual.begin(), qual.end());
    }
    
    string tag_str;
    CuffStrand source_strand = read.source_strand();
    if (source_strand != CUFF_STRAND_UNKNOWN)
    {
        tag_str += source_strand == CUFF_FWD ? "\tXS:A:+" : "\tXS:A:-";
    }
    
    fprintf(sam_out,
            "%s\t%d\t%s\t%d\t255\t%s\t%s\t%d\t0\t%s\t%s%s\n",
            read_name,
            sam_flag,
            ref_name,
            ref_pos + 1,
            cigar_str.c_str(),
            mate_ref_name,
            mate_ref_pos + 1,
            seq.c_str(),
            qual.c_str(),
            tag_str.c_str());
}

void generate_reads(RefSequenceTable& rt,
                    const vector<Scaffold>& ref_mRNAs, 
					AssayProtocol* sequencer,
					int total_frags,
					FILE* sam_frag_out,
					FastqOutfilePair& fastq_out,
					FILE* expr_out,
                    FILE* gtf_out)
{
    RefID last_ref_id = 0;
    int rna_rightmost = 0;
    vector<shared_ptr<ReadHit> > read_chunk; 
	
    if (expr_out)
    {
        fprintf(expr_out, "gene_id\ttranscript_id\tFPKM\trho\tread_cov\tphys_cov\teffective_len\n");
    }
    
	print_sam_header(sam_frag_out, rt);
    
    const FragmentPolicy& frag_policy = sequencer->fragment_policy();
    
	for (size_t i = 0; i < ref_mRNAs.size(); ++i)
	{
        if (last_ref_id &&
            (last_ref_id != ref_mRNAs[i].ref_id() ||
             rna_rightmost <= ref_mRNAs[i].left()))
        {
            sort(read_chunk.begin(), read_chunk.end(), SortReads());
            
            for (size_t j = 0; j < read_chunk.size(); ++j)
            {
                print_aligned_read(*read_chunk[j], rt, sam_frag_out);
            }
            
            read_chunk.clear();
        }
        
		int num_frags_for_mRNA = ref_mRNAs[i].alpha() * total_frags;
        
        vector<bool> covered_by_read(ref_mRNAs[i].length(), false);
		for (size_t j = 0; j < num_frags_for_mRNA; ++j)
		{
			ReadsForFragment reads;
			
			if (sequencer->next_reads(ref_mRNAs[i], reads))
			{
                const ReadHit& left = *(reads.front()); 
                const ReadHit& right = *(reads.back()); 
				read_chunk.push_back(reads.front());
				read_chunk.push_back(reads.back());
                
				print_fastq_pair(reads, fastq_out);
                
                for (size_t k = 0; k < left.read_len(); ++k)
                {
                    covered_by_read[k + left.source_transcript_offset()] = true;
                }
                
                for (size_t k = 0; k < right.read_len(); ++k)
                {
                    covered_by_read[k + right.source_transcript_offset()] = true;
                }
			}
			else
			{
				// TODO: each sequencer should gracefully produce garbage in 
				// cases of bad fragments (i.e. fragment is too small)
			}
		}
        
        int left_cov = -1;
        int transfrag_id = 1;
        for (size_t k = 0; k < covered_by_read.size(); ++k)
        {
            if (covered_by_read[k] && left_cov == -1)
            {
                left_cov = k;
            }
            if (!covered_by_read[k] && left_cov != -1)
            {
                vector<AugmentedCuffOp> transfrag_ops;
                select_genomic_op_range(ref_mRNAs[i].augmented_ops(),
                                        left_cov,
                                        k,
                                        transfrag_ops);
                Scaffold transfrag(ref_mRNAs[i].ref_id(),
                                   ref_mRNAs[i].strand(),
                                   transfrag_ops);
                                   
                string gene_id = ref_mRNAs[i].annotated_gene_id();
                char buf[2048];
                sprintf(buf, "%s_%d", ref_mRNAs[i].annotated_trans_id().c_str(), transfrag_id); 
                vector<string> gtf_recs;
                
                get_scaffold_gtf_records(rt, transfrag, gene_id, buf, gtf_recs);
                foreach (const string& s, gtf_recs)
                {
                    fprintf(gtf_out, "%s\n", s.c_str());
                }
                left_cov = -1;
                transfrag_id++;
            }
        }
        
        if (expr_out)
        {

            double fpkm = 0.0;
            double cov = 0.0;
            double eff_len = ref_mRNAs[i].effective_length(&frag_policy);
            if (eff_len > 0.0)
            {
                fpkm = num_frags_for_mRNA / (eff_len / 1000.0) / (total_frags / 1000000.0);
                cov = (num_frags_for_mRNA * 2 * read_length) / (double)ref_mRNAs[i].effective_length(&frag_policy);
            }

            fprintf(expr_out, 
                    "%s\t%s\t%g\t%g\t%g\t%g\t%g\n", 
                    ref_mRNAs[i].annotated_gene_id().c_str(),
                    ref_mRNAs[i].annotated_trans_id().c_str(),
                    fpkm, 
                    ref_mRNAs[i].rho(),
                    0.0,
                    cov,
                    eff_len);
		}
        
        if (last_ref_id != ref_mRNAs[i].ref_id())
        {
            rna_rightmost = ref_mRNAs[i].right();
        }
        else 
        {
            rna_rightmost = max(rna_rightmost, ref_mRNAs[i].right());
        }
		
		last_ref_id = ref_mRNAs[i].ref_id();
	}
    
    sort(read_chunk.begin(), read_chunk.end(), SortReads());
        
    for (size_t j = 0; j < read_chunk.size(); ++j)
    {
        print_aligned_read(*read_chunk[j], rt, sam_frag_out);
    }
    
    read_chunk.clear();

}

void load_contigs(const string& genome_fasta, 
				  RefSequenceTable& rt, 
				  vector<Scaffold>& source_molecules)
{
	GFastaFile genome(genome_fasta.c_str());
	
	bool last = false;
	do 
	{
		FastaSeq contig;
		genome.getFastaSeq(last, &contig);
		string s = contig.getSeq();
		std::transform(s.begin(), s.end(), s.begin(), (int (*)(int))std::toupper);
		
		vector<AugmentedCuffOp> cuffop(1, AugmentedCuffOp(CUFF_MATCH, 0, s.length()));
		RefID contig_id = rt.get_id(contig.getName(), NULL, s.length());
		source_molecules.push_back(Scaffold(contig_id, CUFF_FWD, cuffop));
		source_molecules.back().seq(s);
	}while(!last);
}

void driver(FILE* sam_out,
			FastqOutfilePair& fastq_out,
			FILE* expr_out,
            FILE* gtf_out,
            FILE* expr_in)
{
	ReadTable it;
	RefSequenceTable rt(true, false);
	
	vector<Scaffold> source_molecules;
	
	if (mrna_gtf != "")
	{        
		FILE* ref_gtf = NULL;
		
		ref_gtf = fopen(mrna_gtf.c_str(), "r");
		if (!ref_gtf)
		{
			fprintf(stderr, "Error: cannot open GTF file %s for reading\n",
					mrna_gtf.c_str());
			exit(1);
		}
		
		load_ref_rnas(ref_gtf, rt, source_molecules);
		if (source_molecules.empty())
		{
			fprintf(stderr, "Error: GTF is empty!\n");
			exit(1);
		}
	}
	else if (genome_fasta != "")
	{
		load_contigs(genome_fasta, rt, source_molecules);
		if (source_molecules.empty())
		{
			fprintf(stderr, "Error: FASTA files is empty!\n");
			exit(1);
		}
	}
	else
	{
		fprintf(stderr, "Error: No source_pool input files defined\n");
		exit(1);
	}
	
	FluxRankAbundancePolicy flux_policy(5e7, -0.6, 9500);
	
    if (expr_in != NULL)
    {
        load_abundances(expr_in, source_molecules);
    }
    else 
    {
        assign_abundances(flux_policy, source_molecules);
    }

	NormalFragments frag_policy(frag_length_mean, 
                                frag_length_std_dev,
                                read_length, frag_length_mean  + 3 * frag_length_std_dev);
	calc_frag_abundances(&frag_policy, source_molecules);
    
	// Set the fragment priming policy, default is uniform random priming.
	if (priming_type == "three_prime")
	{
		shared_ptr<PrimingPolicy> primer = shared_ptr<PrimingPolicy>(new ThreePrimeEndPriming());
		frag_policy.priming_policy(primer);
	}
	
	IlluminaChIPSeqPE seq_policy(read_length, read_length, false);
	
	AssayProtocol* sequencer = new AssayProtocol(frag_policy, seq_policy);
	
	generate_reads(rt, 
                   source_molecules,
				   sequencer,
				   num_fragments,
				   sam_out,
				   fastq_out,
				   expr_out,
                   gtf_out);
	
	delete sequencer;
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
    	
	// seed the random number generator 
	random_seed = time(NULL);
	
	//random_seed = 1111111;
	srand(random_seed);
    
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
	string right_out_fastq_filename = out_prefix + "_2.fq";
	right_fastq_out = fopen(right_out_fastq_filename.c_str(), "w");
	if (!sam_out)
	{
		fprintf(stderr, "Error: cannot open FASTQ file %s for writing\n",
				right_out_fastq_filename.c_str());
		exit(1);
	}
	
	FastqOutfilePair fastq_out(left_fastq_out, right_fastq_out);
    
	FILE* expr_out = NULL;
    FILE* expr_in = NULL;
    
    if (expr_filename != "")
    {
        expr_in = fopen(expr_filename.c_str(), "r");
        if (!expr_in)
        {
            fprintf(stderr, "Error: cannot open expression file %s for reading\n",
                    expr_filename.c_str());
            exit(1);
        }
	}
    else 
    {
        string out_expr_filename = out_prefix + ".simexpr";
        expr_out = fopen(out_expr_filename.c_str(), "w");
        if (!expr_out)
        {
            fprintf(stderr, "Error: cannot open .simexpr file %s for writing\n",
                    out_expr_filename.c_str());
            exit(1);
        }
    }

	

    FILE* frag_gtf_out = NULL;
	string out_gtf_filename = out_prefix + ".transfrags.gtf";
	frag_gtf_out = fopen(out_gtf_filename.c_str(), "w");
	if (!frag_gtf_out)
	{
		fprintf(stderr, "Error: cannot open .simexpr file %s for writing\n",
				out_gtf_filename.c_str());
		exit(1);
	}
    

	driver(sam_out, fastq_out, expr_out, frag_gtf_out, expr_in);

	return 0;
}
