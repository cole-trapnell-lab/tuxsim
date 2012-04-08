/*
 *  scaffolds.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/30/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#include <list>
#include <algorithm>
#include "common.h"
#include "scaffolds.h"
#include "fragments.h"

using namespace std;

bool Scaffold::g_left_lt(const AugmentedCuffOp& lhs,
			 const AugmentedCuffOp& rhs)
{
  return lhs.g_left() < rhs.g_left();
}

bool is_known(const AugmentedCuffOp& op)
{
  return op.opcode != CUFF_UNKNOWN;
}

inline bool has_intron(const Scaffold& scaff)
{
  const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
  for (size_t j = 0; j < ops.size(); ++j)
    {
      if (ops[j].opcode == CUFF_INTRON)
	return true;
    }
  
  return false;
}

// Adds open intervals not covered in the genomic coordinate covered by to_fill
// to the vector gaps.  DOES NOT CLEAR gaps.
void record_gaps(const vector<AugmentedCuffOp>& to_fill, 
		 vector<pair<int, int> >& gaps)
{
  for (size_t i = 1; i < to_fill.size(); ++i)
    {
      if (to_fill[i].g_left() - to_fill[i-1].g_right() > 0)
	{
	  gaps.push_back(make_pair(to_fill[i-1].g_right(), to_fill[i].g_left()));
	}
    }	
}

bool AugmentedCuffOp::compatible(const AugmentedCuffOp& lhs,
				 const AugmentedCuffOp& rhs,
				 bool allow_intron_unknowns)
{
  if (rhs.opcode == CUFF_INTRON)
    {
      if (lhs.opcode == CUFF_INTRON)
	{ 
	  if (lhs != rhs)
	    return false;
	}
      else if (lhs.opcode == CUFF_UNKNOWN)
	{
	  
	}
      else // CUFF_MATCH
	{
	  return false;
	}
    }
  
  if (lhs.opcode == CUFF_INTRON)
    {
      if (rhs.opcode == CUFF_INTRON)
	{
	  if (lhs != rhs)
	    return false;
	}
      else if (lhs.opcode == CUFF_UNKNOWN)
	{
	  
	}
      else // CUFF_MATCH
	{
	  return false;
	}
    }	
  
  return true;
}


bool Scaffold::overlap_in_genome(const Scaffold& lhs, 
				 const Scaffold& rhs, 
				 int overlap_radius)
{
  int ll = lhs.left() - overlap_radius;
  int rr = rhs.right() + overlap_radius;
  int lr = lhs.right() + overlap_radius;
  int rl = rhs.left() - overlap_radius;
  
  if (ll >= rl && ll < rr)
    return true;
  if (lr >rl && lr < rr)
    return true;
  if (rl >= ll && rl < lr)
    return true;
  if (rr > ll && rr < lr)
    return true;
  
  return false;
}

bool intron_op(const AugmentedCuffOp& op)
{
  return op.opcode == CUFF_INTRON;
}

bool Scaffold::strand_agree(const Scaffold& lhs, 
			    const Scaffold& rhs)
{
  bool strand = (lhs.strand() == CUFF_STRAND_UNKNOWN || 
		 rhs.strand() == CUFF_STRAND_UNKNOWN ||
		 lhs.strand() == rhs.strand());
  return strand;
}

bool overlap_in_genome(int ll, int lr, int rl, int rr)
{
  if (ll >= rl && ll < rr)
    return true;
  if (lr > rl && lr < rr)
    return true;
  if (rl >= ll && rl < lr)
    return true;
  if (rr > ll && rr < lr)
    return true;
  return false;
}

int Scaffold::match_length(int left, int right) const
{
  int len = 0;
  size_t curr_op = 0;
  
  while (curr_op != _augmented_ops.size())
    {
      const AugmentedCuffOp& op = _augmented_ops[curr_op];
      
      if (op.opcode == CUFF_MATCH &&
	  ::overlap_in_genome(left, 
			      right, 
			      op.g_left(), 
			      op.g_right()))
	{
	  len += AugmentedCuffOp::match_length(op, left, right);
	}
      if (op.g_left() >= right)
	break;
      ++curr_op;
    }
  
  return len;
}

void Scaffold::clear_hits()
{
  _mates_in_scaff.clear();
}

bool scaff_lt(const Scaffold& lhs, const Scaffold& rhs)
{
  return lhs.left() < rhs.left();
}

bool Scaffold::add_hit(const MateHit* hit)
{
  if (!binary_search(_mates_in_scaff.begin(),
		     _mates_in_scaff.end(),
		     hit))
    {
      _mates_in_scaff.push_back(hit);
    }
  return true;	
}

double Scaffold::effective_length(const FragmentPolicy* frag_policy) const
{
  // Find average of effective lengths
  double eff_len = 0.0;
  double trans_len = length();
  
  // FIXME: replace 1000 with max of frag_policy
  for(int l = 1; l <= trans_len; l++)
    {
      double fp = frag_policy->frag_len_prob(l);
      eff_len += fp * (trans_len - l + 1);
    }
  return eff_len;
}

void Scaffold::insert_true_indels(const vector<AugmentedCuffOp>& indels)
{
  vector<AugmentedCuffOp>::const_iterator low, up;
  low = lower_bound(indels.begin(), indels.end(),
		    AugmentedCuffOp(CUFF_MATCH, _ref_id, _left, 0));
  
  up = upper_bound(indels.begin(), indels.end(),
		   AugmentedCuffOp(CUFF_MATCH, _ref_id, _right, 0));

  vector<AugmentedCuffOp> local_indels;
  if (low == indels.end() || low == up)
    {
      _target_seq = _seq;
      return;
    }

  local_indels.insert(local_indels.begin(), low, up);

  vector<AugmentedCuffOp> temp_ops;
  for (size_t i = 0; i < _augmented_ops.size(); ++i)
    {
      AugmentedCuffOp op = _augmented_ops[i];
      if (op.opcode != CUFF_MATCH)
	{
	  temp_ops.push_back(op);
	  continue;
	}

      bool found_indel = false;
      for (size_t j = 0; j < local_indels.size(); ++j)
	{
	  const AugmentedCuffOp& indel = local_indels[j];
	  if (indel.genomic_offset < op.genomic_offset ||
	      indel.genomic_offset >= op.g_right())
	    continue;

	  if (indel.opcode == CUFF_DEL &&
	      indel.g_right() > op.g_right())
	    continue;

	  found_indel = true;	  
	  AugmentedCuffOp op2 = op;
	  op2.genomic_offset = indel.genomic_offset;
	  if (indel.opcode == CUFF_DEL)
	    op2.genomic_offset += indel.genomic_length;
	  
	  op2.genomic_length = op.g_right() - op2.genomic_offset;
	  
	  op.genomic_length = indel.genomic_offset - op.genomic_offset;
	  if (op.genomic_length > 0)
	    temp_ops.push_back(op);

	  temp_ops.push_back(indel);

	  if (op2.genomic_length > 0)
	    temp_ops.push_back(op2);
	  
	  break;
	}

      if (!found_indel)
	temp_ops.push_back(op);
    }

  _augmented_ops = temp_ops;

  _target_seq = "";
  int start = 0;
  for (size_t i = 0; i < _augmented_ops.size(); ++i)
    {
      const AugmentedCuffOp& op = _augmented_ops[i];
      if (op.opcode == CUFF_MATCH)
	{
	  _target_seq += _seq.substr(start, op.genomic_length);
	  start += op.genomic_length;
	}
      else if (op.opcode == CUFF_INS)
	{
	  for (size_t j = 0; j < op.genomic_length; ++j)
	    _target_seq.push_back('G');
	}
      else if (op.opcode == CUFF_DEL)
	{
	  start += op.genomic_length;
	}
    }

  if (length() != _target_seq.length())
    {
      fprintf(stderr, "length differs! in true indels: %d vs. %d\n",
	      length(), _target_seq.length());
      exit(1);
    }
}

void Scaffold::insert_seq_error_indels(const vector<AugmentedCuffOp>& indels)
{
  string target_seq = "";
  int curr_pos = 0;
  vector<AugmentedCuffOp> temp_ops;
  for (size_t i = 0; i < _augmented_ops.size(); ++i)
    {
      AugmentedCuffOp op = _augmented_ops[i];
      if (op.opcode != CUFF_MATCH)
	{
	  temp_ops.push_back(op);

	  if (op.opcode == CUFF_INS)
	    {
	      target_seq += _target_seq.substr(curr_pos, op.genomic_length);
	      curr_pos += op.genomic_length;
	    }
	  
	  continue;
	}

      bool found_indel = false;
      for (size_t j = 0; j < indels.size(); ++j)
	{
	  AugmentedCuffOp indel = indels[j];
	  if (indel.genomic_offset < curr_pos ||
	      indel.genomic_offset >= curr_pos + op.genomic_length)
	    continue;

	  if (indel.opcode == CUFF_DEL &&
	      indel.genomic_offset + op.genomic_length > curr_pos + op.genomic_length)
	    continue;

	  found_indel = true;
	  AugmentedCuffOp op2 = op;
	  indel.genomic_offset = indel.genomic_offset - curr_pos + op.genomic_offset;
	  op2.genomic_offset = indel.genomic_offset;
	  if (indel.opcode == CUFF_DEL)
	    op2.genomic_offset += indel.genomic_length;
	  
	  op2.genomic_length = op.g_right() - op2.genomic_offset;
	  
	  op.genomic_length = indel.genomic_offset - op.genomic_offset;
	  if (op.genomic_length > 0)
	    {
	      temp_ops.push_back(op);
	      target_seq += _target_seq.substr(curr_pos, op.genomic_length);
	    }


	  temp_ops.push_back(indel);
	  if (indel.opcode == CUFF_INS)
	    {
	      for (size_t k = 0; k < indel.genomic_length; ++k)
		target_seq.push_back('A');
	    }
	  
	  if (op2.genomic_length > 0)
	    {
	      temp_ops.push_back(op2);
	      target_seq += _target_seq.substr(curr_pos + op.genomic_length, op2.genomic_length);
	    }

	  curr_pos += (op.genomic_length + op2.genomic_length);
	  
	  break;
	}

      if (!found_indel)
	{
	  temp_ops.push_back(op);
	  target_seq += _target_seq.substr(curr_pos, op.genomic_length);
	  curr_pos += op.genomic_length;
	}
    }

  _augmented_ops = temp_ops;
  _target_seq = target_seq;

  if (length() != _target_seq.length())
    {
      fprintf(stderr, "length differs! in seq indels: %d vs. %d\n",
	      length(), _target_seq.length());
      exit(1);
    }
}

bool Scaffold::is_in_match(int trans_coord) const
{
  for (size_t i = 0; i < _augmented_ops.size(); ++i)
    {
      const AugmentedCuffOp& op = _augmented_ops[i];
      if (op.opcode == CUFF_MATCH || op.opcode == CUFF_INS)
	{
	  if (trans_coord < op.genomic_length)
	    return op.opcode == CUFF_MATCH;

	  trans_coord -= op.genomic_length;
	}
    }

  assert (0);
  return false;
}

void get_scaffold_gtf_records(const RefSequenceTable& rt, 
                              const Scaffold& scaffold,
                              const string& gene_id,
                              const string& transfrag_id,
                              vector<string>& gtf_recs)
{
  const char* ref_name = rt.get_name(scaffold.ref_id());
  
  const char* strand_str = NULL;
  if (scaffold.strand() == CUFF_STRAND_UNKNOWN)
    strand_str = ".";
  else if (scaffold.strand() == CUFF_FWD)
    strand_str = "+";
  else
    strand_str = "-";
  
  int score = 1000;
  
  char buf[2048];	
  
  sprintf(buf, 
	  "%s\tCufflinks\ttranscript\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\";\n",
	  ref_name,
	  scaffold.left() + 1,
	  scaffold.right(), // GTF intervals are inclusive on both ends, but ours are half-open
	  score,
	  strand_str,
	  gene_id.c_str(),
	  transfrag_id.c_str());
  gtf_recs.push_back(buf);
  
  int exon_num = 1;
  for (size_t op_id = 0; op_id < scaffold.augmented_ops().size(); ++op_id)
    {
      const AugmentedCuffOp& op = scaffold.augmented_ops()[op_id];
      if (op.opcode == CUFF_MATCH || op.opcode == CUFF_UNKNOWN)
	{
	  const char* type = op.opcode == CUFF_MATCH ? "exon" : "missing_data";
	  
	  sprintf(buf, 
		  "%s\tCufflinks\t\%s\t%d\t%d\t%d\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\";\n",
		  ref_name,
		  type,
		  op.g_left() + 1,
		  op.g_right(), // GTF intervals are inclusive on both ends, but ours are half-open
		  score,
		  strand_str,
		  gene_id.c_str(),
		  transfrag_id.c_str());
	  gtf_recs.push_back(buf);
	  
	  exon_num++;
	}
      //gff_recs.push_back(buf);
    }
}
