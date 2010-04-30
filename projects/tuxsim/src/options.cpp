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

}

// If supplied options are incomptabile, prints an error and exits(1)
void validate_options()
{
    bool invalid = false;
    
    if (fastadir == "")
    {
        fprintf(stderr, "Error: input.fastadir must be set\n");
        invalid = true;
    }
    
    // For now, a GTF is required to use TuxSim, which only does RNA-Seq
    if (source_gtf == "")
    {
        fprintf(stderr, "source_pool.source_gtf must be set\n");
        invalid = true;
    }
    
    if (num_fragments <=0)
    {
        fprintf(stderr, "sequencing.num_fragments must be a positive integer\n");
        invalid = true;
    }
    
    if (read_length <= 0)
    {
        fprintf(stderr, "sequencing.num_fragments must be a positive integer\n");
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
        ;
        
        options_description config_file_options("Configuration file options");
        config_file_options.add_options()
            ("input.fasta_dir", value(&fastadir), "")
            ("source_pool.source_gtf", value(&source_gtf), "")
            ("fragment_length.mean", value<double>(&frag_length_mean)->default_value(200), "")
            ("fragment_length.std_dev", value<double>(&frag_length_std_dev)->default_value(40), "")
            //("sequencing.read_type", value(&read_type), "")
            ("sequencing.read_length", value<int>(&read_length)->default_value(76), "")
            ("sequencing.num_fragments", value<int>(&num_fragments)->default_value(20000000), "")
            ("output.prefix", value(&out_prefix), "")
            //("fragment_length.distribution", value(&frag_dist), "")
        ;
        
        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden options");
        hidden.add_options()
        ("input-file", value< vector<string> >(), "input file")
        ;
        
        options_description cmdline_options;
        cmdline_options.add(generic).add(hidden);
        
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
            int cnt = vm.count("input-file");
            vector<string> cfg = vm["input-file"].as< vector<string> >();
            if (cfg.size() != 1)
            {
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
    }
    catch(std::exception& e)
    {
        cerr << e.what() << endl;
    }
    
//    string ofile;
//    string macrofile, libmakfile;
//    bool t_given = false;
//    bool b_given = false;
//    string mainpackage;
//    string depends = "deps_file";
//    string sources = "src_file";
//    string root = ".";
//    
//    options_description desc("Allowed options");
//    desc.add_options()
//    // First parameter describes option name/short name
//    // The second is parameter to option
//    // The third is description
//    ("help,h", "print usage message")
//    ("output,o", value(&ofile), "pathname for output")
//    ("macrofile,m", value(&macrofile), "full pathname of macro.h")
//    ("two,t", bool_switch(&t_given), "preprocess both header and body")
//    ("body,b", bool_switch(&b_given), "preprocess body in the header context")
//    ("libmakfile,l", value(&libmakfile), 
//     "write include makefile for library")
//    ("mainpackage,p", value(&mainpackage), 
//     "output dependency information")
//    ("depends,d", value(&depends), 
//     "write dependencies to <pathname>")
//    ("sources,s", value(&sources), "write source package list to <pathname>")
//    ("root,r", value(&root), "treat <dirname> as project root directory")
//    ;
    return 0;
}