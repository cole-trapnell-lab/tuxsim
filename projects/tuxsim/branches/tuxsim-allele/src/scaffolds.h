#ifndef SCAFFOLDS_H
#define SCAFFOLDS_H
/*
 *  scaffolds.h
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/30/09.
 *  Copyright 2009 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <set>
#include <vector>

#include "common.h"
#include "hits.h"
#include <boost/thread.hpp>

using namespace std;

enum CuffOpCode
{
    CUFF_MATCH,
    CUFF_INTRON,
    CUFF_UNKNOWN,
    CUFF_INS,
    CUFF_DEL,
};

class FragmentPolicy;
class Mismatch;

struct AugmentedCuffOp 
{
    AugmentedCuffOp(const CuffOpCode O, RefID ref_id, int g_off, int g_len) 
	: opcode(O),
    _ref_id(ref_id),
    genomic_offset(g_off),
    genomic_length(g_len) 
    {
        assert (genomic_length >= 0);
    }
    
    RefID ref_id() const { return _ref_id; }
    int g_left() const { return genomic_offset; }
    int g_right() const { return genomic_offset + genomic_length; }
    
    static bool overlap_in_genome(const AugmentedCuffOp& lhs,
                                  const AugmentedCuffOp& rhs)
    {
        if (lhs._ref_id != rhs._ref_id)
            return false;
        
        if (lhs.g_left() >= rhs.g_left() && lhs.g_left() < rhs.g_right())
            return true;
        if (lhs.g_right() > rhs.g_left() && lhs.g_right() < rhs.g_right())
            return true;
        if (rhs.g_left() >= lhs.g_left() && rhs.g_left() < lhs.g_right())
            return true;
        if (rhs.g_right() > lhs.g_left() && rhs.g_right() < lhs.g_right())
            return true;
        return false;
    }
    
    bool contains(const AugmentedCuffOp& other) const
    {
        if (_ref_id != other._ref_id)
            return false;
        
        if (g_left() <= other.g_left() && g_right() >= other.g_right())
            return true;
        return false;
    }
    
    bool properly_contains(const AugmentedCuffOp& other) const
    {
        if (_ref_id != other._ref_id)
            return false;
        
        if ((g_left() < other.g_left() && g_right() >= other.g_right()) ||
            (g_left() <= other.g_left() && g_right() > other.g_right()))
            return true;
        return false;
    }
    
    static int match_length(const AugmentedCuffOp& op, int left, int right)
    {
        int len = 0;
        int left_off = op.g_left();
        if (op.opcode == CUFF_MATCH)
        {
            if (left_off + op.genomic_length > left && left_off < right)
            {
                if (left_off > left)
                {
                    if (left_off + op.genomic_length <= right + 1)
                        len += op.genomic_length;
                    else
                        len += right - left_off;
                }
                else
                {
                    if (left_off + op.genomic_length <= right + 1)
                        len += (left_off + op.genomic_length - left);
                    else
                        return right - left;
                }
            }
        }
        return len;
    }
    
    static bool compatible(const AugmentedCuffOp& lhs,
                           const AugmentedCuffOp& rhs,
                           bool allow_intron_unknowns);
    
    bool operator==(const AugmentedCuffOp& rhs) const
    {
		bool result;
		//allele
         if(opcode == rhs.opcode &&
			_ref_id == rhs._ref_id &&
			genomic_offset == rhs.genomic_offset &&
			genomic_length == rhs.genomic_length)
			 result = true;
		 else
			 result = false;
		 return result;
	}
    
    bool operator<(const AugmentedCuffOp& rhs) const
    {
        if (_ref_id != rhs._ref_id)
            return _ref_id < rhs._ref_id;
        
        if (genomic_offset != rhs.genomic_offset)
            return genomic_offset < rhs.genomic_offset;
        
        if (genomic_length != rhs.genomic_length)
            return genomic_length > rhs.genomic_length;
        
        if (opcode != rhs.opcode)
            return opcode == CUFF_MATCH;
        return false;
    }
    
    bool operator!=(const AugmentedCuffOp& rhs) const
    {
        return !(*this == rhs);
    }
    
    CuffOpCode opcode;
    RefID _ref_id;
    int genomic_offset;
    int genomic_length;
};

class Scaffold
{
    void cuff_ops_from_cigar(vector<AugmentedCuffOp>& ops,
                             const vector<CigarOp>& cig,
                             int& g_left)
	{
        for (size_t i = 0; i < cig.size(); ++i)
	    {
            assert(cig[i].length >= 0);
            switch(cig[i].opcode)
            {
                case MATCH:
                    ops.push_back(AugmentedCuffOp(CUFF_MATCH, _ref_id, g_left, cig[i].length));
                    g_left += cig[i].length;
                    break;
                case REF_SKIP:
                    ops.push_back(AugmentedCuffOp(CUFF_INTRON, _ref_id, g_left, cig[i].length));
                    g_left += cig[i].length;
                    break;
                case SOFT_CLIP:
                    g_left += cig[i].length;
                    break;
                default:
                    assert(false);
                    break;
            }
	    }
	}
    
    RefID _ref_id;
    
public:
    Scaffold() :
    _ref_id(0), 
    _strand(CUFF_STRAND_UNKNOWN), 
    _classcode(0),	  
    _rho(0.0),
    _alpha(0.0),
    _frags_generated(0) {}
    
    // For manually constructing scaffolds, for example when a reference is 
    // available
    Scaffold(RefID ref_id, CuffStrand strand, const vector<AugmentedCuffOp>& ops)
    : _ref_id(ref_id), 
    _augmented_ops(ops), 
    _strand(strand),
    _classcode(0),
    _rho(0.0),
    _alpha(0.0),
    _frags_generated(0)
    {
        _right = _augmented_ops.back().g_right();
        _left = _augmented_ops.front().g_left();
        _has_intron = has_intron(*this);
        
        assert(!_has_intron || _strand != CUFF_STRAND_UNKNOWN);
    }
    
    //int get_id() const { return _id; }
    
    int left() const { return _left; }
    
    int right() const  { return _right; }
    
    const string& annotated_trans_id() const { return _annotated_trans_id; }
    void annotated_trans_id(const string& ann_name) { _annotated_trans_id = ann_name; }
    
    const string& annotated_gene_id() const { return _annotated_gene_id; }
    void annotated_gene_id(const string& ann_name) { _annotated_gene_id = ann_name; }
    
    const string& annotated_gene_name() const { return _annotated_gene_name; }
    void annotated_gene_name(const string& ann_name) { _annotated_gene_name = ann_name; }
    
    const string& annotated_protein_id() const { return _annotated_protein_id; }
    void annotated_protein_id(const string& ann_name) { _annotated_protein_id = ann_name; }	
    
    const string& annotated_tss_id() const { return _annotated_tss_id; }
    void annotated_tss_id(const string& ann_name) { _annotated_tss_id = ann_name; }	
    
    const string& nearest_ref_id() const { return _nearest_ref_id; }
    void nearest_ref_id(const string& ann_name) { _nearest_ref_id = ann_name; }	
    
    char nearest_ref_classcode() const { return _classcode; }
    void nearest_ref_classcode(char cc) { _classcode = cc; }
    
    bool has_intron() const { return _has_intron; }
    
    CuffStrand strand() const { return _strand; }
    void strand(CuffStrand strand) { _strand = strand; }
    
    // Could we merge lhs and rhs together?
    static bool compatible(const Scaffold& lhs, 
                           const Scaffold& rhs, 
                           bool allow_intron_unknowns);
    
    static bool strand_agree(const Scaffold& lhs, 
                             const Scaffold& rhs);
    
    
    const vector<const MateHit*>& mate_hits() const { return _mates_in_scaff; }
    RefID ref_id() const { return _ref_id; }
    
    bool strictly_contains(const Scaffold& other) const;
    
    bool contains(const Scaffold& other) const
    {
        return (left() <= other.left() && right() >= other.right());
    }
    
    int match_length(int left, int right) const;
    
    double effective_length(const FragmentPolicy* frag_policy) const;
    
    int length() const
    {
        int len = 0;
        
        for (size_t j = 0; j < _augmented_ops.size(); ++j)
        {
            if (_augmented_ops[j].opcode == CUFF_MATCH ||
                _augmented_ops[j].opcode == CUFF_INS)
                len += _augmented_ops[j].genomic_length;
        }
        
        return len;
    }
    
    static bool g_left_lt(const AugmentedCuffOp& lhs,
                          const AugmentedCuffOp& rhs);
	
    const vector<AugmentedCuffOp>& augmented_ops() const { return _augmented_ops; }
    
  
  double rho() const { return _rho; }
  void rho(double r) { _rho = r; }
  
  double alpha() const { return _alpha; }
  void alpha(double a) { _alpha = a; }

  void insert_true_mismatches(const vector<Mismatch>& mismatches);
  void insert_seq_error_mismatches(const vector<Mismatch>& mismatches);
  void insert_true_indels(const vector<AugmentedCuffOp>& indels);
  void insert_seq_error_indels(const vector<AugmentedCuffOp>& indels);
  
  bool is_in_match(int trans_coord) const;

  int gseq_id() const { return _gseq_id; }
  void gseq_id(int id) { _gseq_id = id; }

  // NOTE: sequences are always stored in forward strand space, they are NOT 
  // reverse complemented
  const string& seq() const { return _seq; }
  void seq(const string& s) { _seq = s; }
  
  const string& target_seq() const
  {
    if (_target_seq.empty())
      return _seq;
    
    return _target_seq;
  }
  
	//allele
  string target_seq()
  {
    if (_target_seq.empty())
      return _seq;
    
    return _target_seq;
  }
  
  int num_frags_generated() const { return _frags_generated; }
  void num_frags_generated(int f) { _frags_generated = f; }
    
 private: 
  
  static bool has_intron(const Scaffold& scaff)
  {
    const vector<AugmentedCuffOp>& ops = scaff.augmented_ops();
    for (size_t j = 0; j < ops.size(); ++j)
      {
	if (ops[j].opcode == CUFF_INTRON)
	  return true;
      }

    return false;
  }

static bool overlap_in_genome(const Scaffold& lhs, 
			      const Scaffold& rhs, 
			      int overlap_radius);

vector<pair< int, int> > gaps() const
    {
        vector<pair<int,int> > g;
        const vector<AugmentedCuffOp>& ops = augmented_ops();
        for (size_t j = 0; j < ops.size(); ++j)
        {
            if (ops[j].opcode == CUFF_INTRON)
            {
                g.push_back(make_pair(ops[j].g_left(), ops[j].g_right()));
            }
        }
        return g;
    }
    
    inline bool has_unknown() const
    {
        const vector<AugmentedCuffOp>& ops = augmented_ops();
        for (size_t j = 0; j < ops.size(); ++j)
        {
            if (ops[j].opcode == CUFF_UNKNOWN)
                return true;
        }
        return false;
    }
    
    void clear_hits();
    bool add_hit(const MateHit*);
    
private: 
    typedef vector<AugmentedCuffOp> OpList;
    
    vector<const MateHit*> _mates_in_scaff;
    
    int _left;
    int _right;
    bool _has_intron; 
    string _seq;
    string _target_seq;
    
    vector<AugmentedCuffOp> _augmented_ops;
    CuffStrand _strand;
    
    string _annotated_trans_id;
    string _annotated_gene_id;
    string _annotated_gene_name;
    string _annotated_protein_id;
    string _annotated_tss_id;
    string _nearest_ref_id;
    char _classcode;
    
    double _rho;
    double _alpha;
    int _frags_generated;

  int _gseq_id;
};

bool scaff_lt(const Scaffold& lhs, const Scaffold& rhs);

bool overlap_in_genome(int ll, int lr, int rl, int rr);

struct StructurallyEqualScaffolds
{
    bool operator()(const Scaffold& lhs, const Scaffold& rhs)
    {
        return lhs.ref_id() == rhs.ref_id() && 
        lhs.augmented_ops() == rhs.augmented_ops();
    }
};

void get_scaffold_gtf_records(const RefSequenceTable& rt, 
                              const Scaffold& scaffold,
                              const string& gene_id,
                              const string& transfrag_id,
                              vector<string>& gtf_recs);
#endif
