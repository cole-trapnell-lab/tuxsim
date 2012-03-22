/*
 *  sequencing.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "sequencing.h"

using namespace std;

bool IlluminaChIPSeqPE::reads_for_fragment(const LibraryFragment& frag, 
                                           ReadsForFragment& reads)
{
  assert (frag.source_seq);
  const Scaffold& mRNA = *(frag.source_seq);
  int frag_length = frag.end - frag.start;;
  if (frag_length < _left_len || frag_length < _right_len)
    return false;
  
  _next_fragment_id++;
  
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
  
  const string& target_seq = frag.source_seq->target_seq();
  string left_seq, right_seq;
  
  left_seq = target_seq.substr(frag.start, _left_len);
  right_seq = target_seq.substr(frag.end - _right_len, _right_len);
  
  bool reverse_strand_frag = _bool_generator();
  
  shared_ptr<ReadHit> left_read(new ReadHit());
  shared_ptr<ReadHit> right_read(new ReadHit());
  
  CuffStrand strand = CUFF_STRAND_UNKNOWN;
  
  *left_read = ReadHit(frag.source_seq->ref_id(),
		       _next_fragment_id,
		       left_read_ops.front().g_left(),
		       left_read_cigar,
		       true,
		       strand,
		       frag.source_seq->ref_id(),
		       right_read_ops.front().g_left(),
		       0,
		       0,
		       0, //mRNA.annotated_trans_id()
		       frag.start);
  
  *right_read = ReadHit(frag.source_seq->ref_id(),
			_next_fragment_id,
			right_read_ops.front().g_left(),
			right_read_cigar,
			true,
			strand,
			frag.source_seq->ref_id(),
			left_read_ops.front().g_left(),
			0,
			0,
			0, //mRNA.annotated_trans_id()
			frag.end - _right_len);

  // daehwan - for debug purposes
  if (left_read->read_len() != _left_len ||
      right_read->read_len() != _right_len)
    {
      fprintf(stderr, "this should not happen!!\n");
      exit(1);
    }

  
  if (_strand_specific || left_read->has_intron() || right_read->has_intron())
    {
      left_read->source_strand(frag.source_seq->strand());
      right_read->source_strand(frag.source_seq->strand());
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
		  // daehwan - for debug purposes
		  fprintf(stderr, "%d) trans(%d-%d), genomic(%d-%d)\n",
			  ++count, start, end, genomic_left, genomic_right);
		  assert(false);
		  exit(1);
		  return false;
		}
	      
	      genomic_left = curr_op.genomic_offset + (start - rna_counter);
            }
	  if (genomic_right == -1 &&
	      end <= rna_counter + curr_op.genomic_length)
            {
	      if (curr_op.opcode == CUFF_INS)
		{
		  // daehwan - for debug purposes
		  fprintf(stderr, "%d) trans(%d-%d), genomic(%d-%d)\n",
			  ++count, start, end, genomic_left, genomic_right);
		  assert(0);
		  exit(1);
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
	    out_ops.push_back(curr_op);
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
