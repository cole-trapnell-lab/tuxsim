/*
 *  fragments.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "fragments.h"

void NormalFragments::next_fragment(const Scaffold& molecule,
                                    LibraryFragment& fragment)
{
    int frag_length;
    do {
        frag_length = _length_generator();
    }while (frag_length < _min_frag_length || frag_length > _max_frag_length);

    int frag_start = (molecule.length() - frag_length) * _uniform_generator();
    frag_start = max(frag_start, 0);
    
    assert (frag_start < molecule.length());
    
    int frag_end = frag_start + frag_length;
    frag_end = min(frag_end, molecule.length());
    //assert (frag_end <= molecule.length());
    
    fragment.source_seq = &molecule;
    fragment.start = frag_start;
    fragment.end = frag_end;
}