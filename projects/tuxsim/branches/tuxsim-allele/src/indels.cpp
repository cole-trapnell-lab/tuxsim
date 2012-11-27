#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "indels.h"

void generate_true_indels(const vector<AugmentedCuffOp>& exons,
			  vector<AugmentedCuffOp>& indels)
{
  if (indel_true_diff_per_bases <= 0)
    return;
  
  for (size_t i = 0; i < exons.size(); ++i)
    {
      const AugmentedCuffOp& exon = exons[i];
      
      int random_number = rand() % indel_true_diff_per_bases;
      if (random_number < exon.genomic_length)
	{
	  CuffOpCode opcode;
	  if (rand() % 2 == 0)
	    opcode = CUFF_INS;
	  else
	    opcode = CUFF_DEL;
	  
	  int length = rand() % 3 + 1;
	  int pos = exon.g_left() + rand() % exon.genomic_length;
	  if (opcode == CUFF_DEL &&
	      pos + length >= exon.g_right())
	    continue;
	  
	  indels.push_back(AugmentedCuffOp(opcode, exon._ref_id, pos, length));
	}
    }
}

void print_indels(FILE* indels_out,
		  RefSequenceTable& rt,
		  const vector<AugmentedCuffOp>& indels)
{
  for (size_t i = 0; i < indels.size(); ++i)
    {
      const AugmentedCuffOp& indel = indels[i];
      const char* ref_name = rt.get_name(indel.ref_id());
      if (indel.opcode == CUFF_INS)
	{
	  fprintf(indels_out, "I %s %d %d\n",
		  ref_name, indel.genomic_offset, indel.genomic_length);
	}
      else
	{
	  fprintf(indels_out, "D %s %d %d\n",
		  ref_name, indel.genomic_offset, indel.genomic_length);
	}
    }
}
