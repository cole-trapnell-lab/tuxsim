#ifndef FRAGMENTS_H
#define FRAGMENTS_H
/*
 *  fragments.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/random/normal_distribution.hpp>

#include "scaffolds.h"

struct LibraryFragment
{
	const Scaffold* source_seq;
	
	// Fragment specified as an open interval [start, end)
	int start;
	int end;	
};

class FragmentPolicy
{
public:
	virtual bool next_fragment(const Scaffold& molecule,
							   LibraryFragment& fragment) = 0;
};

class NormalFragments : public FragmentPolicy
{
	
	// This is a typedef for a random number generator.
	// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
	typedef boost::minstd_rand base_generator_type;
	
	typedef variate_generator<base_generator_type&, boost::uniform_01<> > uniform_generator_type;
	typedef variate_generator<base_generator_type&, normal_distribution<> > normal_generator_type;
	
public:
	
	NormalFragments(int mean_frag_length, 
                    int frag_length_sd, 
                    int min_frag_length,
                    int max_frag_length) :
        _base_generator(base_generator_type(random_seed)),
        _uniform_generator(uniform_generator_type(_base_generator, uniform_01<>())),
        _length_generator(normal_generator_type(_base_generator, 
                                                normal_distribution<>(mean_frag_length, frag_length_sd))),
        _min_frag_length(min_frag_length),
        _max_frag_length(max_frag_length)
	{
	}
	
	virtual bool next_fragment(const Scaffold& molecule,
							   LibraryFragment& fragment);
	
private:
	base_generator_type _base_generator;
	uniform_generator_type _uniform_generator;
	normal_generator_type _length_generator;
    
    int _min_frag_length;
    int _max_frag_length;
};
#endif
