#ifndef MISMATCHES_H
#define MISMATCHES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <map>
#include <vector>
#include <cassert>

using namespace std;

#include "common.h"
#include "scaffolds.h"
#include "hits.h"

struct Mismatch
{
Mismatch(RefID ref_id, int g_off) :
  _ref_id(ref_id),
    _genomic_offset(g_off)
  {
  }

  RefID ref_id() const { return _ref_id; }
  int genomic_offset() const { return _genomic_offset; }

  bool operator==(const Mismatch& rhs) const
  {
    return (_ref_id == rhs._ref_id &&
	    _genomic_offset == rhs._genomic_offset);
  }
  
  bool operator<(const Mismatch& rhs) const
  {
    if (_ref_id != rhs._ref_id)
      return _ref_id < rhs._ref_id;
    
    if (_genomic_offset != rhs._genomic_offset)
      return _genomic_offset < rhs._genomic_offset;

    return false;
  }
  
  bool operator!=(const Mismatch& rhs) const
  {
    return !(*this == rhs);
  }

    
  RefID _ref_id;
  int _genomic_offset;
};

void generate_true_mismatches(const vector<AugmentedCuffOp>& exons,
			      vector<Mismatch>& mismatches);

void print_mismatches(FILE* mismatches_out,
		      RefSequenceTable& rt,
		      const vector<Mismatch>& mismatches);


#endif
