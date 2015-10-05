/*
 *  options.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/14/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

//#include <getopt.h>
#include <string>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

#include "options.h"

#pragma GCC visibility push(default) 

using namespace boost;
using namespace boost::program_options;
using namespace std;

void print_usage()
{
    fprintf(stderr, "tuxsim <sim.cfg>\n");
    fprintf(stderr, "--------------------------------------------\n");
}

// If supplied options are incomptabile, prints an error and exits(1)
void validate_options()
{
    bool invalid = false;

	//allele
	if(allele_simulator){
		if (mrna_gtf == "" || fastadir == ""){
			fprintf(stderr, "Error: for allele_simulator both mrna_gtf and fastadir must be set\n");
			invalid = true;
		}
		if(vcf_table == ""){
			fprintf(stderr, "Error: input.vcf_table must be set\n");
			invalid = true;
		}			
	}
	
    if (mrna_gtf != "" && fastadir == "")
    {
        fprintf(stderr, "Error: input.fastadir must be set\n");
        invalid = true;
    }
    
    //    // For now, a GTF is required to use TuxSim, which only does RNA-Seq
    //    if (mrna_gtf == "")
    //    {
    //        fprintf(stderr, "source_pool.mrna_gtf must be set\n");
    //        invalid = true;
    //    }
    
	if (priming_type == "uniform_random")
	{
		
	}
	else if (priming_type == "three_prime")
	{
		
	}
	else
	{
		fprintf(stderr, "Error: fragment.priming must be one of [uniform_random,three_prime]\n");
		exit(1);
	}
	
    if (num_fragments <=0)
    {
		fprintf(stderr, "Error: sequencing.num_fragments must be a positive integer\n");
        invalid = true;
    }
    
    if (read_length <= 0)
    {
        fprintf(stderr, "Error: sequencing.num_fragments must be a positive integer\n");
        invalid = true;
    }
    
    if (invalid)
        exit(1);
}

int parse_options(int argc, char** argv)
{
    try
    {
        options_description generic("Command line options");
        generic.add_options()
        ("help,h", "print usage message")
        //("output.prefix", value(&out_prefix), "")
        ("expression,e", value(&expr_filename), "Load mRNA expression values from input file instead of generating them")
        ;
        
        options_description config_file_options("Configuration file options");
        config_file_options.add_options()
	  ("input.fasta_dir", value(&fastadir), "")
			//allele
  	  ("input.allele_simulator", value<bool>(&allele_simulator)->default_value(false), "")
	  ("input.sort_by_position", value<bool>(&sort_by_position)->default_value(true), "")
	  ("input.only_phased_reads", value<bool>(&only_phased_reads)->default_value(false), "")
	  ("output.prefix", value(&out_prefix), "")
	  ("source_pool.mrna_gtf", value(&mrna_gtf)->default_value(""), "")
	  ("source_pool.genome_fasta", value(&genome_fasta)->default_value(""), "")
	  //allele
	  ("source_pool.vcf_table", value(&vcf_table)->default_value(""), "")
	  ("fragment.priming", value(&priming_type)->default_value("uniform_random"), "")
	  ("fragment.length.mean", value<double>(&frag_length_mean)->default_value(200), "")
	  ("fragment.length.std_dev", value<double>(&frag_length_std_dev)->default_value(40), "")
	  //("sequencing.read_type", value(&read_type), "")
	  ("sequencing.read_length", value<int>(&read_length)->default_value(75), "")
	  ("sequencing.num_fragments", value<int>(&num_fragments)->default_value(20000000), "")
	  ("sequencing.max_edit_dist", value<int>(&max_edit_dist)->default_value(0), "")
      ("sequencing.single_end", value<bool>(&single_end)->default_value(false), "")	
	  ("mismatch.true_diff_per_bases", value<int>(&mismatch_true_diff_per_bases)->default_value(0), "")
	  ("mismatch.seq_error_per_bases", value<int>(&mismatch_seq_error_per_bases)->default_value(0), "")
	  ("indel.true_diff_per_bases", value<int>(&indel_true_diff_per_bases)->default_value(0), "")
	  ("indel.seq_error_per_bases", value<int>(&indel_seq_error_per_bases)->default_value(0), "")
          ("ase.silenced_allele", value<string>(&silenced_allele)->default_value(""), "")
          ("ase.silenced_fraction", value<double>(&silenced_fraction)->default_value(1.0), "")
          ("ase.allele_proportion_natural_stdev", value<double>(&allele_proportion_natural_stdev)->default_value(0.1), "")
	  //("fragment_length.distribution", value(&frag_dist), "")
        ;
        
        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden options");
        hidden.add_options()
        ("input-file", value< vector<string> >(), "input file")
        ;
        
        options_description cmdline_options;
        cmdline_options.add(generic).add(config_file_options).add(hidden);
        
        positional_options_description p;
        p.add("input-file", -1);
        
        variables_map vm;
        store(command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);
        
        if (vm.count("help")) {  
            cout << cmdline_options << "\n";
            return 0;
        }
        
        if (vm.count("input-file"))
        {
            vm.count("input-file");
            vector<string> cfg = vm["input-file"].as< vector<string> >();
            if (cfg.size() != 1)
            {
                print_usage();
                return 1;
            }
            string cfg_file = cfg[0];
            
            ifstream ifs(cfg_file.c_str());
            
            if (!ifs.good())
            {
				fprintf(stderr, "Error: cannot load config file %s\n", cfg_file.c_str());
                return 1;
            }
            
            store(parse_config_file(ifs, config_file_options), vm);
            notify(vm);
            
            validate_options();
        }
        else
        {
            print_usage();
            return 1;
        }
    }
    catch(std::exception& e)
    {
        cerr << e.what() << endl;
    }
    
    return 0;
}
