/*
 *  fragments.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "fragments.h"


bool UniformRandomPriming::next_priming_site(const Scaffold& molecule,
											 int& priming_site)
{
	priming_site = molecule.length() * _uniform_generator();
	
	return true;
}


bool ThreePrimeEndPriming::next_priming_site(const Scaffold& molecule,
											 int& priming_site)
{
	priming_site =  molecule.length();
	return true;
}

bool NormalFragments::next_fragment(const Scaffold& molecule,
                                    LibraryFragment& fragment)
{
    int frag_length;
	
	if (_min_frag_length >= molecule.length())
	{
		return false;
	}
	
	int frag_end;
	
	bool valid_priming_site;
	
	do {
		valid_priming_site = _priming_policy->next_priming_site(molecule, 
																frag_end);
	}
	while (valid_priming_site && frag_end < _min_frag_length);
	
	if (!valid_priming_site)
	{
		return false;
	}
	
    do {
        frag_length = _length_generator();
		if (frag_end - frag_length < 0)
			frag_length = frag_end;
		
    }while (frag_length < _min_frag_length || frag_length > _max_frag_length);

    int frag_start = frag_end - frag_length;
    frag_start = max(frag_start, 0);
    
    assert (frag_start < molecule.length());
    
	//assert (frag_end <= molecule.length());
    
    fragment.source_seq = &molecule;
    fragment.start = frag_start;
    fragment.end = frag_end;
	
	return true;
}