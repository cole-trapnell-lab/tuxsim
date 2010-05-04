/*
 *  common.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 3/31/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/foreach.hpp>
#define foreach         BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH

using namespace std;
using namespace boost;

extern string priming_type;
extern double frag_length_mean;
extern double frag_length_std_dev;
extern int num_fragments;
extern int read_length;
extern int num_reads;
extern string out_prefix;

extern string genome_fasta;

// Both fastadir and mrna_gtf must be set, pending resolution of #169
extern string fastadir;
extern string mrna_gtf;

extern int random_seed;

int parseInt(int lower, 
			 const char *errmsg, 
			 void (*print_usage)());

float parseFloat(float lower, 
				 float upper, 
				 const char *errmsg, 
				 void (*print_usage)());

void reverse_complement(string& seq);


// From SAMtools (Li et al):
/*! the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! the mate is unmapped */
#define BAM_FMUNMAP        8
/*! the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! this is read1 */
#define BAM_FREAD1        64
/*! this is read2 */
#define BAM_FREAD2       128
/*! not primary alignment */
#define BAM_FSECONDARY   256
/*! QC failure */
#define BAM_FQCFAIL      512
/*! optical or PCR duplicate */
#define BAM_FDUP        1024