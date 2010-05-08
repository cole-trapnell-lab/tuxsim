/*
 *  abundances.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "abundances.h"

void assign_expression_ranks(const vector<Scaffold>& ref_mRNAs,
							 vector<unsigned int>& expr_rank)
{
	for (size_t i = 0; i < expr_rank.size(); ++i)
	{
		expr_rank[i] = (unsigned int)i;
	}
	
	random_shuffle(expr_rank.begin(), expr_rank.end());
}

// FIXME: should use effective length
void calc_frag_abundances(vector<Scaffold>& source_molecules)
{
	double fragment_pool_size = 0.0;
	
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		fragment_pool_size += source_molecules[i].length() * source_molecules[i].rho();
	}
	
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		double frag_n = source_molecules[i].length() * source_molecules[i].rho();
		source_molecules[i].alpha(frag_n / fragment_pool_size);
	}
}