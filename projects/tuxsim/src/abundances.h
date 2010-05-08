#ifndef ABUNDANCES_H
#define ABUNDANCES_H

/*
 *  abundances.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <vector>
#include "common.h"
#include "scaffolds.h"

void assign_expression_ranks(const vector<Scaffold>& ref_mRNAs,
							 vector<unsigned int>& expr_rank);

// Assigns the transcriptome wide relative abundances
template<class RankedBasedAbundancePolicy>
void assign_abundances(const RankedBasedAbundancePolicy& rank_based,
					   vector<Scaffold>& source_molecules)
{
	vector<unsigned int> expr_rank(source_molecules.size(), 0);
	vector<double> expr_rho(source_molecules.size(), 0.0);
	
	assert(expr_rho.size() == expr_rank.size());
	
	assign_expression_ranks(source_molecules, expr_rank);
	
	double total_mols = 0.0;
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		int rank = expr_rank[i];
		double rho = rank_based.rho(rank);
		expr_rho[i] = rho;
		total_mols += rho;
	}
	
	assert (total_mols > 0.0);
	
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		expr_rho[i] /= total_mols;
		source_molecules[i].rho(expr_rho[i]);
	}
}

void calc_frag_abundances(vector<Scaffold>& source_molecules);


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
#endif
