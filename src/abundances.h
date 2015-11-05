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
#include <map>
#include <ctime>

#include "common.h"
#include "fragments.h"
#include "scaffolds.h"

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

void assign_expression_ranks(vector<unsigned int>& expr_rank);

// Assigns the transcriptome wide relative abundances
template<class RankedBasedAbundancePolicy>
void assign_abundances(const RankedBasedAbundancePolicy& rank_based,
					   vector<Scaffold>& source_molecules)
{
        int n_ranks = source_molecules.size();
        if (allele_simulator) n_ranks /= 2;
	vector<unsigned int> expr_rank(n_ranks, 0);
	vector<double> expr_rho(source_molecules.size(), 0.0);

	assign_expression_ranks(expr_rank);
	
	double total_mols = 0.0;

        if (!allele_simulator) {
            for (int i = 0; i < (int) source_molecules.size(); ++i) {
                int rank = expr_rank[i];
                double rho = rank_based.rho(rank);
                expr_rho[i] = rho;
                total_mols += rho;
            }
        } else {
            int j = 0;
            map<string, double> base_rho;
            for (int i = 0; i < (int) source_molecules.size(); ++i) {
                string transcript_id = source_molecules[i].annotated_trans_id();
                string base_transcript = transcript_id.substr(0, transcript_id.size()-2);
                string allele = transcript_id.substr(transcript_id.size()-1, 1);
                if (allele == "P") {
                    int rank = expr_rank[j];
                    double rho = rank_based.rho(rank);
                    base_rho[base_transcript] = rho;
                    if (allele == silenced_allele)
                        rho *= (1.0 - silenced_fraction);

                    expr_rho[i] = rho;
                    total_mols += rho;
                    j++;
                }
            }

            boost::minstd_rand rng;
            rng.seed(time(NULL));
            boost::normal_distribution<> allele_prop_dist(0.5, allele_proportion_natural_stdev);
            boost::variate_generator<boost::minstd_rand&, boost::normal_distribution<> > sample_allele_proportion(rng, allele_prop_dist);
            

            for (int i = 0; i < (int) source_molecules.size(); ++i) {
                string transcript_id = source_molecules[i].annotated_trans_id();
                string base_transcript = transcript_id.substr(0, transcript_id.size()-2);
                string allele = transcript_id.substr(transcript_id.size()-1, 1);
                if (allele == "M") {
                    double rho = base_rho[base_transcript];

                    double allele_proportion = sample_allele_proportion();
                    rho = max(rho * 2.0 * allele_proportion, 0.0);
                    if (allele == silenced_allele)
                        rho *= (1.0 - silenced_fraction);

                    expr_rho[i] = rho;
                    total_mols += rho;
                }
            }
        }

	assert (total_mols > 0.0);
	
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		expr_rho[i] /= total_mols;
		source_molecules[i].rho(expr_rho[i]);
	}
}

void load_abundances(FILE* expr_in, vector<Scaffold>& source_molecules);

void calc_frag_abundances(const FragmentPolicy* frag_policy,
                          vector<Scaffold>& source_molecules);


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
