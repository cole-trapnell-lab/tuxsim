#ifndef INDELS_H
#define INDELS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

#include "common.h"
#include "scaffolds.h"
#include "hits.h"

void generate_true_indels(const vector<AugmentedCuffOp>& exons,
			  vector<AugmentedCuffOp>& indels);

void print_indels(FILE* indels_out,
		  RefSequenceTable& rt,
		  const vector<AugmentedCuffOp>& indels);

#endif
