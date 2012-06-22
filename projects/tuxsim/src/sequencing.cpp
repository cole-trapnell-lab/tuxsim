/*
 *  sequencing.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "common.h"
#include "sequencing.h"

#include "GStr.h"
#include "GFaSeqGet.h"
#include "GFastaFile.h"
#include "gtf_tracking.h"

#include "mismatches.h"

using namespace std;

int bowtie_sam_extra(int gseq_id, const ReadHit& rh, GFastaHandler& gfasta, vector<string>& fields);

bool IlluminaChIPSeqPE::reads_for_fragment(const LibraryFragment& frag, 
                                           ReadsForFragment& reads,
					   GFastaHandler& gfasta)
{
  assert (frag.source_seq);
  Scaffold mRNA = *(frag.source_seq);

  int frag_length = frag.end - frag.start;;
  if (frag_length < _left_len || frag_length < _right_len)
    return false;

  if (mismatch_seq_error_per_bases > 0)
    {
      vector<Mismatch> mismatches;
      int random_number = rand() % mismatch_seq_error_per_bases;
      if (random_number < frag_length)
	{
	  mismatches.push_back(Mismatch(mRNA.ref_id(), random_number));
	}
      
      mRNA.insert_seq_error_mismatches(mismatches);
    }  
  
  if (indel_seq_error_per_bases > 0)
    {
      vector<AugmentedCuffOp> indels;
      int random_number = rand() % indel_seq_error_per_bases;
      if (random_number < frag_length)
	{
	  bool pass = true;
	  static const int edge = 5;
	  if (random_number < edge)
	    {
	      random_number = edge;
	      pass = false;
	    }
	  else if (random_number >= frag_length - edge)
	    {
	      random_number = frag_length - edge - 1;
	      pass = false;
	    }
	  else if (random_number < _left_len &&
		   _left_len - random_number < edge)
	    {
	      random_number = edge * 2;
	      pass = false;
	    }
	  
	  if (random_number > frag_length - _right_len - edge &&
		   random_number - (frag_length - _right_len) < edge)
	    pass = false;

	  if (pass)
	    {
	      CuffOpCode opcode;
	      if (rand() % 2 == 0)
		opcode = CUFF_INS;
	      else
		opcode = CUFF_DEL;

	      int length = rand() % 3 + 1;
	      indels.push_back(AugmentedCuffOp(opcode, mRNA.ref_id(), frag.start + random_number, length));
	      
	      mRNA.insert_seq_error_indels(indels);
	    }
	}
    }
    
    vector<AugmentedCuffOp> frag_ops;
    const vector<AugmentedCuffOp>& rna_ops = mRNA.augmented_ops();
    
    vector<AugmentedCuffOp> left_read_ops;
    bool result = select_genomic_op_range(rna_ops, frag.start, frag.start + _left_len, left_read_ops);
    if (!result)
        return false;
    
    assert (!left_read_ops.empty());
    
    vector<CigarOp> left_read_cigar;
    cuff_op_to_cigar(left_read_ops, left_read_cigar);
    assert (!left_read_cigar.empty());
    
    vector<AugmentedCuffOp> right_read_ops;
    result = select_genomic_op_range(rna_ops, frag.end - _right_len, frag.end, right_read_ops);
    if (!result)
        return false;
    
    assert (!right_read_ops.empty());
    
    vector<CigarOp> right_read_cigar;
    cuff_op_to_cigar(right_read_ops, right_read_cigar);
    assert (!right_read_cigar.empty());
    
    const string& target_seq = mRNA.target_seq();
    string left_seq, right_seq;
    
    left_seq = target_seq.substr(frag.start, _left_len);
    right_seq = target_seq.substr(frag.end - _right_len, _right_len);
    
    bool reverse_strand_frag = _bool_generator();
    
    shared_ptr<ReadHit> left_read(new ReadHit());
    shared_ptr<ReadHit> right_read(new ReadHit());
    
    CuffStrand strand = CUFF_STRAND_UNKNOWN;
    
    _next_fragment_id++;
    
    *left_read = ReadHit(mRNA.ref_id(),
                         _next_fragment_id,
                         left_read_ops.front().g_left(),
                         left_read_cigar,
                         true,
                         strand,
                         mRNA.ref_id(),
                         right_read_ops.front().g_left(),
                         0,
                         0,
                         0, //mRNA.annotated_trans_id()
                         frag.start);
    
    *right_read = ReadHit(mRNA.ref_id(),
                          _next_fragment_id,
                          right_read_ops.front().g_left(),
                          right_read_cigar,
                          true,
                          strand,
                          mRNA.ref_id(),
                          left_read_ops.front().g_left(),
                          0,
                          0,
                          0, //mRNA.annotated_trans_id()
                          frag.end - _right_len);
    
    if (left_read->read_len() != _left_len ||
        right_read->read_len() != _right_len)
    {
        fprintf(stderr, "this should not happen!!\n");
        exit(1);
    }
    
    if (_strand_specific || left_read->has_intron() || right_read->has_intron())
    {
        left_read->source_strand(mRNA.strand());
        right_read->source_strand(mRNA.strand());
    }
    
    int base_flag = BAM_FPAIRED | BAM_FPROPER_PAIR;
    int left_flag = base_flag;
    int right_flag = base_flag;
    
    string right_seq_rc = right_seq;            
    reverse_complement(right_seq_rc);
    right_flag |= BAM_FREVERSE;
    left_flag |= BAM_FMREVERSE;
    
    left_read->seq(left_seq);
    right_read->seq(right_seq_rc);
    
    if (reverse_strand_frag)
    {
        left_flag |= BAM_FREAD2;
        right_flag |= BAM_FREAD1;
    }
    else
    {
        left_flag |= BAM_FREAD1;
        right_flag |= BAM_FREAD2;
    }
  
  reads.push_back(left_read);
  reads.push_back(right_read);
  
  left_read->sam_flag(left_flag);
  right_read->sam_flag(right_flag);
  
  char buf[2048];
  sprintf(buf,"%llu",left_read->insert_id());
  
  left_read->name(buf);
  right_read->name(buf);
  
  left_read->qual(string(left_read->seq().length(), 'I'));
  right_read->qual(string(right_read->seq().length(), 'I'));

  vector<string> left_aux_fields;
  int left_edit_dist = bowtie_sam_extra(mRNA.gseq_id(), *left_read, gfasta, left_aux_fields);
  left_read->aux_sam_fields(left_aux_fields);

  vector<string> right_aux_fields;
  int right_edit_dist = bowtie_sam_extra(mRNA.gseq_id(), *right_read, gfasta, right_aux_fields);
  right_read->aux_sam_fields(right_aux_fields);

  if (left_edit_dist > max_edit_dist || right_edit_dist > max_edit_dist)
    return false;

  return true;
}

bool select_genomic_op_range(const vector<AugmentedCuffOp>& src_ops,
                             int start,
                             int end,
                             vector<AugmentedCuffOp>& out_ops)
{
    out_ops.clear();
    if (src_ops.empty())
        return false;
    
    int genomic_left = -1;
    int genomic_right = -1;
    int rna_counter = 0;
    
    // daehwan - for debug purposes
    static int count = 0;
#if 0
    fprintf(stderr, "\n\ntranscript\n");
    int rna_coord = 0;
    for (size_t i = 0; i < src_ops.size(); ++i)
    {
        const AugmentedCuffOp& op = src_ops[i];
        fprintf(stderr, "\topcode(%d): %d-%d (%d-%d in transcriptome)\n",
                op.opcode, op.g_left(), op.g_right(),
                rna_coord, rna_coord + op.genomic_length);
        
        if (op.opcode == CUFF_MATCH || op.opcode == CUFF_INS)
            rna_coord += op.genomic_length;
    }
    
    fprintf(stderr, "\nrna coord: %d-%d\n\n",
            start, end);
#endif
    
    size_t curr_op_idx = 0;
    while (curr_op_idx < src_ops.size())
    {
        // the calculation below is tracked in RNA coordinates, but the result
        // is saved under the projection to genomic coordinates.
        const AugmentedCuffOp& curr_op = src_ops[curr_op_idx];
        if (curr_op.opcode == CUFF_MATCH || curr_op.opcode == CUFF_INS)
        {
            if (genomic_left == -1 &&
                start < rna_counter + curr_op.genomic_length)
            {
                if (curr_op.opcode == CUFF_INS)
                {
                    /*
                     fprintf(stderr, "%d) trans(%d-%d), genomic(%d-%d)\n",
                     ++count, start, end, genomic_left, genomic_right);
                     */
                    
                    return false;
                }
                
                genomic_left = curr_op.genomic_offset + (start - rna_counter);
            }
            if (genomic_right == -1 &&
                end <= rna_counter + curr_op.genomic_length)
            {
                if (curr_op.opcode == CUFF_INS)
                {
                    /*
                     fprintf(stderr, "%d) trans(%d-%d), genomic(%d-%d)\n",
                     ++count, start, end, genomic_left, genomic_right);
                     */
                    
                    return false;
                }
                
                genomic_right = curr_op.genomic_offset + (end - rna_counter);
                break;
            }
            
            rna_counter += curr_op.genomic_length;
        }
        
        curr_op_idx++;
    }
    
    int curr_genomic_off = src_ops.front().g_left();
    curr_op_idx = 0;
    
#if 0
    fprintf(stderr, "\tgenome coord: %d-%d\n",
            genomic_left, genomic_right);
#endif
    
    while (curr_op_idx < src_ops.size() &&
           curr_genomic_off < src_ops.back().g_right() && 
           curr_genomic_off < genomic_right)
    {
        const AugmentedCuffOp& curr_op = src_ops[curr_op_idx++];
        
        if (curr_genomic_off + curr_op.genomic_length > genomic_left &&
            curr_genomic_off < genomic_right)
        {
            if (curr_op.opcode != CUFF_INS || genomic_left < curr_genomic_off)
            {
                out_ops.push_back(curr_op);
            }
        }
        
        if (curr_op.opcode != CUFF_INS)
            curr_genomic_off += curr_op.genomic_length;
    }
    
    if (out_ops.empty())
        return false;
    
    if (out_ops.front().opcode != CUFF_MATCH || out_ops.back().opcode != CUFF_MATCH)
    {
        fprintf(stderr, "%d) front or back is not MATCH\n", ++count);
        exit(1);
        assert (0);
        return false;
    }
    
    if (genomic_right <= out_ops.back().g_left() ||
        genomic_right > out_ops.back().g_right())
    {
        fprintf(stderr, "%d) right is out of range\n", ++count);
        exit(1);
        assert(0);
        return false;
    }
    
    int right_clip = out_ops.back().g_right() - genomic_right;
    out_ops.back().genomic_length -= right_clip;
    
    if (genomic_left < out_ops.front().g_left() ||
        genomic_left >= out_ops.front().g_right())
    {
        fprintf(stderr, "%d) left is out of range\n", ++count);
        exit(1);
        assert(0);
        return false;
    }
    
    int left_clip = genomic_left - out_ops.front().g_left();
    out_ops.front().genomic_length -= left_clip;
    out_ops.front().genomic_offset = genomic_left;
    
    return true;
}

void cuff_op_to_cigar(const vector<AugmentedCuffOp>& cuff_ops,
                      vector<CigarOp>& cigar)
{
    cigar.clear();
    for (size_t i = 0; i < cuff_ops.size(); ++i)
    {
        switch(cuff_ops[i].opcode)
        {
            case CUFF_MATCH:
                cigar.push_back(CigarOp(MATCH, cuff_ops[i].genomic_length));
                break;
                
            case CUFF_INTRON:
                cigar.push_back(CigarOp(REF_SKIP, cuff_ops[i].genomic_length));
                break;
                
            case CUFF_INS:
                cigar.push_back(CigarOp(INS, cuff_ops[i].genomic_length));
                break;
                
            case CUFF_DEL:
                cigar.push_back(CigarOp(DEL, cuff_ops[i].genomic_length));
                break;
                
            default:
                assert(false);
                break;
        }
    }
}

void str_appendInt(string& str, int64_t v)
{
  char int_str[32] = {0};
  sprintf(int_str, "%ld", v);
  str += int_str;
}

int bowtie_sam_extra(int gseq_id, const ReadHit& rh, GFastaHandler& gfasta, vector<string>& fields)
{
  static GFaSeqGet* buf_faseq = NULL;
  static int buf_gseq_id = -1;

  GFaSeqGet* faseq = NULL;
  if (gseq_id != buf_gseq_id)
    {
      buf_gseq_id = gseq_id;
      faseq = buf_faseq = gfasta.fetch(gseq_id);
    }
  else
    {
      faseq = buf_faseq;
    }

  if (!faseq)
    return -1;
  
  int length = 1;
  const char* ref_str = faseq->subseq(0, length);

  if (!ref_str)
    return -1;

  size_t pos_seq = 0;
  size_t pos_mismatch = 0;
  size_t pos_ref = rh.left();
  size_t mismatch = 0;
  size_t N_mismatch = 0;
  size_t num_gap_opens = 0;
  size_t num_gap_conts = 0;

  static const int bowtie2_min_score = -10;
  static const int bowtie2_max_penalty = 6;
  static const int bowtie2_min_penalty = 2;
  static const int bowtie2_penalty_for_N = 1;
  static const int bowtie2_read_gap_open = 5;
  static const int bowtie2_read_gap_cont = 3;
  static const int bowtie2_ref_gap_open = 5;
  static const int bowtie2_ref_gap_cont = 3;

  int AS_score = 0;
  
  const vector<CigarOp>& cigars = rh.cigar();
  string seq = rh.seq();
  if (rh.sam_flag() & BAM_FREVERSE)
    reverse_complement(seq);

  const string& qual = rh.qual();
  string AS = "AS:i:";
  string MD = "MD:Z:";

  for (size_t i = 0; i < cigars.size(); ++i)
    {
      CigarOp cigar = cigars[i];
      switch(cigar.opcode)
	{
	case MATCH:
	  {
	    const char* ref_seq = ref_str + pos_ref + 1;
	    pos_ref += cigar.length;
	    
	    for (size_t j = 0; j < cigar.length; ++j)
	      {
		char ref_nt = ref_seq[j];
		if (seq[pos_seq] != ref_nt)
		  {
		    ++mismatch;

		    if (pos_seq < qual.length())
		      {
			if (seq[pos_seq] == 'N' || ref_nt == 'N')
			  {
			    AS_score -= (int)bowtie2_penalty_for_N;
			  }
			else
			  {
			    float penalty = bowtie2_min_penalty + (bowtie2_max_penalty - bowtie2_min_penalty) * min((int)(qual[pos_seq] - '!'), 40) / 40.0;
			    AS_score -= (int)penalty;
			  }
		      }

		    str_appendInt(MD, (int)pos_mismatch);
		    MD.push_back((char)ref_nt);
		    pos_mismatch = 0;
		  }
		else
		  {
		    if (ref_nt == 'N')
		      {
			++N_mismatch;
			AS_score -= (int)bowtie2_penalty_for_N;
		      }

		    ++pos_mismatch;
		  }

		++pos_seq;
	      }
	  }
	  break;

	case INS:
	  {
	    pos_seq += cigar.length;

	    AS_score -= bowtie2_read_gap_open;
	    AS_score -= (int)(bowtie2_read_gap_cont * cigar.length);

	    num_gap_opens += 1;
	    num_gap_conts += cigar.length;
	  }
	  break;
	  
	case DEL:
	  {
	    AS_score -= bowtie2_ref_gap_open;
	    AS_score -= (int)(bowtie2_ref_gap_cont * cigar.length);

	    num_gap_opens += 1;
	    num_gap_conts += cigar.length;
	      
	    const char* ref_seq = ref_str + pos_ref + 1;
	    pos_ref += cigar.length;

	    str_appendInt(MD, (int)pos_mismatch);
	    MD.push_back('^');
	    for (size_t k = 0; k < cigar.length; ++k)
	      MD.push_back((char)ref_seq[k]);
	    
	    pos_mismatch = 0;
	  }
	  break;

	case REF_SKIP:
	  {
	    pos_ref += cigar.length;
	  }
	  break;

	default:
	  break;
	}
    }

  str_appendInt(AS, AS_score);
  fields.push_back(AS);

  string XM = "XM:i:";
  str_appendInt(XM, (int)mismatch);
  fields.push_back(XM);

  string XO = "XO:i:";
  str_appendInt(XO, (int)num_gap_opens);
  fields.push_back(XO);

  string XG = "XG:i:";
  str_appendInt(XG, (int)num_gap_conts);
  fields.push_back(XG);

  str_appendInt(MD, (int)pos_mismatch);
  fields.push_back(MD);

  string NM = "NM:i:";
  int edit_dist = mismatch + num_gap_conts;
  str_appendInt(NM, (int)edit_dist);
  fields.push_back(NM);

  return edit_dist;
}
