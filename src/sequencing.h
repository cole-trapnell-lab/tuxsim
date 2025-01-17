#ifndef SEQUENCING_H
#define SEQUENCING_H
/*
 *  sequencing.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/12/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "fragments.h"
#include "hits.h"

typedef vector<boost::shared_ptr<ReadHit> > ReadsForFragment;

class GFastaHandler;

class SequencingPolicy
{
public:
  virtual ~SequencingPolicy() {}
  virtual bool reads_for_fragment(const LibraryFragment& frag, 
				  ReadsForFragment& reads,
				  GFastaHandler& gfasta) = 0;
};

/*******************************************************************************
 Sequencing Policies
 *******************************************************************************/
struct IlluminaChIPSeqPE : public SequencingPolicy
{
  
  typedef boost::minstd_rand base_generator_type;
  typedef variate_generator<base_generator_type&,
    boost::uniform_smallint<> > bool_generator_type;
  
 IlluminaChIPSeqPE(int left_read_len, 
		   int right_read_len,
		   bool strand_specific) :
  _left_len(left_read_len),
    _right_len(right_read_len),
    _strand_specific(strand_specific),
    _next_fragment_id(0),
    _base_generator(base_generator_type(random_seed)),
    _bool_generator(bool_generator_type(_base_generator, uniform_smallint<>(0,1)))
      {}
  
  bool reads_for_fragment(const LibraryFragment& frag, 
			  ReadsForFragment& reads,
			  GFastaHandler& gfasta);
 private:
  int _left_len;
  int _right_len;
  bool _strand_specific;
  InsertID _next_fragment_id;
  base_generator_type _base_generator;
  bool_generator_type _bool_generator;
};

// Selects a list of (clipped) cuff_ops from an RNA, given an interval
// in RNA coordinates	
bool select_genomic_op_range(const vector<AugmentedCuffOp>& src_ops,
                             int start,
                             int end,
                             vector<AugmentedCuffOp>& out_ops);

void cuff_op_to_cigar(const vector<AugmentedCuffOp>& cuff_ops,
                      vector<CigarOp>& cigar);
#endif
