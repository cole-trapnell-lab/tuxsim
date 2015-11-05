/*
 *  common.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 3/31/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "common.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>

string priming_type = "uniform_random";
double frag_length_mean = 200;
double frag_length_std_dev = 40;
int num_fragments = 20000000;
int read_length = 75;
int max_edit_dist = 0;
bool single_end = false;
string out_prefix;

string genome_fasta;

int mismatch_true_diff_per_bases = 0;
int mismatch_seq_error_per_bases = 0;

int indel_true_diff_per_bases = 0;
int indel_seq_error_per_bases = 0;

string fastadir;
string mrna_gtf;
//allele
string vcf_table;
bool allele_simulator = false;
bool sort_by_position = true;
bool only_phased_reads = false;
string expr_filename;
string silenced_allele;
double silenced_fraction = 1.0;
double allele_proportion_natural_stdev = 0.1;

int random_seed;

using namespace std;

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */

int parseInt(int lower, const char *errmsg, void (*print_usage)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            cerr << errmsg << endl;
            print_usage();
            exit(1);
        }
        return (int32_t)l;
    }
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */
float parseFloat(float lower, float upper, const char *errmsg, void (*print_usage)()) {
    float l;
    l = (float)atof(optarg);
	
    if (l < lower) {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    if (l > upper)
    {
        cerr << errmsg << endl;
        print_usage();
        exit(1);
    }
	
    return l;
	
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}

void str_appendInt(string& str, int64_t v)
{
    char int_str[32] = {0};
    sprintf(int_str, "%ld", v);
    str += int_str;
}

// This could be faster.
void reverse_complement(string& seq)
{
    //fprintf(stderr,"fwd: %s\n", seq.c_str());
    for (string::size_type i = 0; i < seq.length(); ++i)
    {
        switch(seq[i])
        {
            case 'A' : seq[i] = 'T'; break;
            case 'T' : seq[i] = 'A'; break;
            case 'C' : seq[i] = 'G'; break;
            case 'G' : seq[i] = 'C'; break;
            default: 
            {
                seq[i]   = 'N'; break;
            }
        }
    }
    reverse(seq.begin(), seq.end());
    //fprintf(stderr, "rev: %s\n", seq.c_str());
}

//allele
void splitString(const string& str,vector<string>& subStrs,const string& delimiter)
{
	// Skip delimiter at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiter,0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiter,lastPos);
	subStrs.clear();
	
	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a subStr, add it to the vector.
		subStrs.push_back(str.substr(lastPos,pos - lastPos));
		// Skip delimiter.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiter,pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiter,lastPos);
	}
}