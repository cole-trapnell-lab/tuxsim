/*
 *  options.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/14/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

const char *short_options = "s:m:";

static struct option long_options[] = {
	{"inner-dist-mean",			required_argument,       0,          'm'},
	{"inner-dist-stddev",		required_argument,       0,          's'},
	{0, 0, 0, 0} // terminator
};

void print_usage()
{
	fprintf(stderr, "tuxsim v%s\n", PACKAGE_VERSION); 
	fprintf(stderr, "-----------------------------\n"); 
	
	//NOTE: SPACES ONLY, bozo
    fprintf(stderr, "Usage:   tuxsim [options] <ref_transcripts.gtf> <fasta_di> <out_sample_name>\n");
	fprintf(stderr, "Options:\n\n");
	fprintf(stderr, "-m/--frag-length-mean         the mean fragment length                         [ default:     200 ]\n");
	fprintf(stderr, "-s/--frag-length-std-dev      the standard deviation of the fragment length    [ default:     50 ]\n");
}

int parse_options(int argc, char** argv)
{
    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
        switch (next_option) {
			case -1:     /* Done with options. */
				break;
			case 'm':
				frag_length_mean = (uint32_t)parseInt(-200, "-m/--frag-length-mean arg must be at least -200", print_usage);
				break;
			case 's':
				frag_length_std_dev = (uint32_t)parseInt(0, "-s/--frag-length-std-dev arg must be at least 0", print_usage);
				break;
			default:
				print_usage();
				return 1;
        }
    } while(next_option != -1);
	
	
	return 0;
}