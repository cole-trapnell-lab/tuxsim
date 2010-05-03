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

double frag_length_mean = 200;
double frag_length_std_dev = 40;
int num_fragments = 20000000;
int read_length = 76;
string out_prefix;

string fastadir;
string source_gtf;

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
                int a = 5;
                seq[i]   = 'N'; break;
            }
        }
    }
    reverse(seq.begin(), seq.end());
    //fprintf(stderr, "rev: %s\n", seq.c_str());
}