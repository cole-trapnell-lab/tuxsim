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
		const AugmentedCuffOp&  op = _augmented_ops[curr_op];
		
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
