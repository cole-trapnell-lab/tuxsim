/*
 *  abundances.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <stdlib.h>

#include "abundances.h"

void load_abundances(FILE* expr_in, vector<Scaffold>& source_molecules)
{
	//vector<unsigned int> expr_rank(source_molecules.size(), 0);
	vector<double> expr_rho(source_molecules.size(), 0.0);
    
    char orig_buf[2048];
    map<string, double> rhos;
    
    double total_rho = 0.0;
    while (fgets(orig_buf, 2048, expr_in))
	{
        const char* buf = orig_buf;
        const char* _gene_name = strsep((char**)&buf,"\t");
            
        if (!_gene_name)
            return;
        char gene_name[2048];
        strncpy(gene_name, _gene_name, 2047); 
        
        if (!strcmp(gene_name, "gene_id"))
            continue;
        
        const char* rna_name = strsep((char**)&buf,"\t");
        if (!rna_name)
            return;
        
        const char* FPKM = strsep((char**)&buf,"\t");
        if (!FPKM)
            return;
        
        const char* _rho = strsep((char**)&buf,"\t");
        if (!_rho)
            return;
        
        const char* read_cov =  strsep((char**)&buf,"\t");
        if (!read_cov)
            return;
        
        const char* phys_cov = strsep((char**)&buf,"\t");
        if (!phys_cov)
            return;
        
        const char* tss_id = strsep((char**)&buf,"\t");
        if (!tss_id)
            return;

        double rho = strtod(_rho, NULL); 
        rhos[rna_name] = rho;
        
        total_rho += rho;
    }
    
    if (total_rho == 0.0)
        return;
    
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
        const string& name = source_molecules[i].annotated_trans_id();
        map<string, double>::const_iterator ci = rhos.find(name);
        if (ci == rhos.end() || total_rho == 0.0)
        {
            source_molecules[i].rho(0.0);
        }
        else 
        {
            double rho = ci->second / total_rho;
            if (rho < 0.0)
            {
                fprintf(stderr, "Error: rho must be >= 0.0!\n");
                exit(1);
            }
            source_molecules[i].rho(rho);
        }
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

void calc_frag_abundances(const FragmentPolicy* frag_policy,
                          vector<Scaffold>& source_molecules)
{
	double fragment_pool_size = 0.0;
	
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		fragment_pool_size += source_molecules[i].effective_length(frag_policy) * source_molecules[i].rho();
	}
	
    if (fragment_pool_size == 0.0)
        return;
    
	for (size_t i = 0; i < source_molecules.size(); ++i)
	{
		double frag_n = source_molecules[i].effective_length(frag_policy) * source_molecules[i].rho();
        double new_alpha = frag_n / fragment_pool_size;
        if (new_alpha > 1.0 || new_alpha < 0.0 || isnan(new_alpha) || isinf(new_alpha))
        {
            fprintf(stderr, "Error: bad alpha (%lg) for RNA!\n", new_alpha);
            fprintf(stderr, "\t rho : %lg, effective_length : %lg for RNA!\n", source_molecules[i].rho(), source_molecules[i].effective_length(frag_policy));
            exit(1);
        }
        
		source_molecules[i].alpha(new_alpha);
	}
}