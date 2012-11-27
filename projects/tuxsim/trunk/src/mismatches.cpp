#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "mismatches.h"

void generate_true_mismatches(const vector<AugmentedCuffOp>& exons,
			      vector<Mismatch>& mismatches)
{
  if (mismatch_true_diff_per_bases <= 0)
    return;
  
  for (size_t i = 0; i < exons.size(); ++i)
    {
      const AugmentedCuffOp& exon = exons[i];
      
      int random_number = rand() % mismatch_true_diff_per_bases;
      if (random_number < exon.genomic_length)
	{
	  int pos = exon.g_left() + rand() % exon.genomic_length;
	  mismatches.push_back(Mismatch(exon.ref_id(), pos));
	}
    }
}

void print_mismatches(FILE* mismatches_out,
		      RefSequenceTable& rt,
		      const vector<Mismatch>& mismatches)
{
  for (size_t i = 0; i < mismatches.size(); ++i)
    {
      const Mismatch& mismatch = mismatches[i];
      const char* ref_name = rt.get_name(mismatch.ref_id());
      fprintf(mismatches_out, "%s %d\n",
	      ref_name, mismatch.genomic_offset());
    }
}
