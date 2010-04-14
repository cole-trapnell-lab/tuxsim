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

using namespace std;
using namespace boost;

extern int frag_length_mean;
extern int frag_length_std_dev;

extern int max_intron_length;

extern const char* fastadir;

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