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

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
typedef boost::minstd_rand base_generator_type;

struct LibraryFragment
{
	const Scaffold* source_seq;
	
	// Fragment specified as an open interval [start, end)
	int start;
	int end;	
};

class PrimingPolicy
{
public:
	virtual bool next_priming_site(const Scaffold& molecule,
								   int& priming_site) = 0;
};

class UniformRandomPriming : public PrimingPolicy
{
	
	typedef variate_generator<base_generator_type&, boost::uniform_01<> > uniform_generator_type;
	
public:
	UniformRandomPriming() : 
		_base_generator(base_generator_type(random_seed)),
		_uniform_generator(uniform_generator_type(_base_generator, uniform_01<>()))
	{
	}
	bool next_priming_site(const Scaffold& molecule,
						   int& priming_site);
	
private:
	base_generator_type _base_generator;
	uniform_generator_type _uniform_generator;
};

class ThreePrimeEndPriming : public PrimingPolicy
{
	bool next_priming_site(const Scaffold& molecule,
						   int& priming_site);
};

class FragmentPolicy
{
public:
	virtual bool next_fragment(const Scaffold& molecule,
							   LibraryFragment& fragment) = 0;
};

class NormalFragments : public FragmentPolicy
{
	typedef variate_generator<base_generator_type&, normal_distribution<> > normal_generator_type;
	
public:
	
	NormalFragments(int mean_frag_length, 
                    int frag_length_sd, 
                    int min_frag_length,
                    int max_frag_length) :
        _base_generator(base_generator_type(random_seed)),
        _length_generator(normal_generator_type(_base_generator, 
                                                normal_distribution<>(mean_frag_length, frag_length_sd))),
        _min_frag_length(min_frag_length),
        _max_frag_length(max_frag_length)
	{
		_priming_policy = shared_ptr<PrimingPolicy>(new UniformRandomPriming());
	}
	
	virtual bool next_fragment(const Scaffold& molecule,
							   LibraryFragment& fragment);
	
	void priming_policy(shared_ptr<PrimingPolicy> policy ) 
	{ 
		_priming_policy = policy; 
	}
	
private:
	base_generator_type _base_generator;
	normal_generator_type _length_generator;
	
    shared_ptr<PrimingPolicy> _priming_policy;
    int _min_frag_length;
    int _max_frag_length;
};
#endif
