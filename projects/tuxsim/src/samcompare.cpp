/*
 *  samcmp.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/30/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

//#include <getopt.h>
#include <string>

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iterator>

#include <list>
#include <set>
#include <map>

#include "common.h"

#include "hits.h"

using namespace boost;
using namespace boost::program_options;
using namespace std;

string ref_sam_name;
string target_sam_name;

void print_usage()
{
    fprintf(stderr, "samcompare  [options] <ref.sam> <target.sam>\n");
    fprintf(stderr, "--------------------------------------------\n");
}

int parse_options(int argc, char** argv)
{
    try
    {
        options_description generic("Command line options");
        generic.add_options()
        ("help,h", "print usage message")
        ;
         
        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden options");
        hidden.add_options()
        ("input-file", value< vector<string> >(), "input file")
        ;
        
        options_description cmdline_options;
        cmdline_options.add(generic).add(hidden);
        
        positional_options_description p;
        p.add("input-file", -1);
        
        variables_map vm;
        store(command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
        notify(vm);
        
        if (vm.count("help")) {  
            cout << cmdline_options << "\n";
            return 0;
        }
        
        if (vm.count("input-file"))
        {
            vector<string> args = vm["input-file"].as< vector<string> >();
            if (args.size() != 2)
            {
                print_usage();
                cerr << generic << endl;
                exit(1);
            }
            ref_sam_name = args[0];
            target_sam_name = args[1];
        }
        else
        {
            print_usage();
            cerr << generic << endl;
            exit(1);
        }
    }
    catch(std::exception& e)
    {
        cerr << e.what() << endl;
    }

    return 0;
}


/*******************************************************************************
 HitBundle is a set of MateHit objects that, were you to look at the interval
 graph of their spanning intervals in genomic coordinates, you'd see a single
 connected component. Note that bundles do not correspond to single transcripts,
 or even single genes
 *******************************************************************************/
class HitBundle
{
public:
	HitBundle() : _leftmost(INT_MAX), _rightmost(-1), _final(false), _id(++_next_id) {}
	~HitBundle() 
	{
	}
	
	int left()   const { return _leftmost;  }
	int right()  const { return _rightmost; }
	int length() const { return _rightmost - _leftmost; }
	
	// Returns true if the hit was added successfully.
	bool add_hit(const MateHit& hit);
	
	const std::vector<MateHit>& hits() const { return _hits; } 
	const std::vector<MateHit>& non_redundant_hits() const { return _non_redundant; } 
	const std::vector<int>& collapse_counts() const { return _collapse_counts; } 
	
	RefID ref_id()  const
	{
		if (!_hits.empty())
			return _hits.front().ref_id();
		else
			return 0;
	}
	
	int id() const { return _id; }
	
	// Adds a Bowtie hit to the open hits buffer.  The Bundle will handle turning
	// the Bowtie hit into a properly mated Cufflinks hit record
	void add_open_hit(shared_ptr<ReadHit> bh);
	
	// Commits any mates still open as singleton hits
	void finalize_open_mates();
	
	
private:
	int _leftmost;
	int _rightmost;
	std::vector<MateHit> _hits;
	std::vector<MateHit> _non_redundant;
	std::vector<int> _collapse_counts;
	bool _final;
	int _id;
	
	static int _next_id;
	
	typedef map<int, list<MateHit> > OpenMates;
	OpenMates _open_mates;
};

int HitBundle::_next_id = 0;

bool HitBundle::add_hit(const MateHit& hit)
{
	if (_final)
		return false;
	
	// Update the bounds on the span
	if (hit.left() < _leftmost)
		_leftmost = hit.left();
	if (hit.right() > _rightmost)
		_rightmost = hit.right();
	
	_hits.push_back(hit);
	return true;
}

void HitBundle::add_open_hit(shared_ptr<ReadHit> bh)
{
	if (bh->partner_ref_id() == 0)
	{
		// This is a singleton, so just make a closed MateHit and
		// continue
		MateHit m(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0);
		add_hit(m);
	}
	else
	{
		OpenMates::iterator mi = _open_mates.find(bh->left());
		
		// Does this read hit close an open mate?
		if (mi == _open_mates.end())
		{
			// No, so add it to the list of open mates, unless we would
			// already have seen it's partner
			if(bh->left() < bh->partner_pos())
			{
				MateHit open_hit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0);
				
				pair<OpenMates::iterator, bool> ret;
				ret = _open_mates.insert(make_pair(bh->partner_pos(), 
                                                   list<MateHit>()));
				
				ret.first->second.push_back(open_hit);
			}
			else
			{
				add_hit(MateHit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0));
			}
		}
		else
		{
			bool found_partner = false;
			// Maybe, see if we can find an ID match in the list of
			// open mates expecting a partner at this position
			for (list<MateHit>::iterator pi = mi->second.begin();
				 pi != mi->second.end();
				 ++pi)
			{
				MateHit& pm = *pi;
				
				if (pm.insert_id() == bh->insert_id())
				{
                    pm.right_alignment(bh);
                    add_hit(pm);
                    mi->second.erase(pi);
                    if (mi->second.empty())
                        _open_mates.erase(mi);
                    
                    found_partner = true;
                    break;
				}
			}
			
			if (!found_partner)
			{
				// If we got here, couldn't actually close any mates with
				// this read hit, so open a new one, unless we can never
				// close this one
				if(bh->left() < bh->partner_pos())
				{
					MateHit open_hit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0);
					
					pair<OpenMates::iterator, bool> ret;
					ret = _open_mates.insert(make_pair(bh->partner_pos(), 
                                                       list<MateHit>()));
					
					ret.first->second.push_back(open_hit);
				}
				else
				{
					add_hit(MateHit(bh->ref_id(), bh, shared_ptr<ReadHit const>(), 0, 0));
				}
			}
		}
	}
}

void HitBundle::finalize_open_mates()
{
	for (OpenMates::iterator itr = _open_mates.begin(); 
		 itr != _open_mates.end(); 
		 ++itr)
	{
		for (list<MateHit>::iterator mi = itr->second.begin(); mi != itr->second.end(); ++mi)
		{
			add_hit(*mi);
		}
	}
    
    sort(_hits.begin(), _hits.end());
}


// Returns groups of all alignments that start at a given position.  For ease
// of comparing two SAM files without blowing up memory.
class PositionBundleFactory
{
public:
	PositionBundleFactory(SAMHitFactory& fac, FILE* hfile)
    : sam_hit_fac(fac), hit_file(hfile), _next_line_num(0) {}
	
	bool next_bundle(HitBundle& bundle_out);
	
	SAMHitFactory& hit_factory() { return sam_hit_fac; } 
	
	void reset() { rewind(hit_file); _next_line_num = 0; }
	
private:
	SAMHitFactory sam_hit_fac;
	FILE* hit_file;
	int _next_line_num;
};

bool PositionBundleFactory::next_bundle(HitBundle& bundle_out)
{
	HitBundle bundle = bundle_out;
	
	if (feof(hit_file))
	{
		return false;
	}
	char bwt_buf[2048];
	
	RefID last_ref_id_seen = 0;
	int last_pos_seen = 0;
	
	off_t curr_pos = ftello(hit_file);
	
	while (fgets(bwt_buf, 2048, hit_file))
	{
		// Chomp the newline
		char* nl = strrchr(bwt_buf, '\n');
		if (nl) *nl = 0;
		
		_next_line_num++;
		
		shared_ptr<ReadHit> bh(new ReadHit());
		
		bool more_hits = sam_hit_fac.get_hit_from_buf(_next_line_num, bwt_buf, *bh, false);
        if (!more_hits)
		{
			continue;
		}
        
        bh->orig_hit_str(bwt_buf);
		
		if (bh->ref_id() == 84696373) // corresponds to SAM "*" under FNV hash. unaligned read record  
			continue;
    
        if (!last_ref_id_seen || 
            (last_ref_id_seen == bh->ref_id() && last_pos_seen == bh->left()))
        {
            bundle.add_open_hit(bh);
        }
        else if (last_ref_id_seen == bh->ref_id() && (bh->left() < last_pos_seen))
        {
            fprintf(stderr, "Error: this SAM file doesn't appear to be correctly sorted!\n");
            fprintf(stderr, "\tcurrent hit is at %s:%d, last one was at %s:%d\n", 
                    sam_hit_fac.ref_table().get_name(bh->ref_id()),
                    bh->left(),
                    sam_hit_fac.ref_table().get_name(last_ref_id_seen),
                    last_pos_seen);
            
            exit(1);
        }
        else
        {
            fseeko(hit_file, curr_pos, SEEK_SET);
            break;
        }
		
		last_ref_id_seen = bh->ref_id();
		last_pos_seen = bh->left();
		
		curr_pos = ftello(hit_file);
	}
	
	bundle.finalize_open_mates();
	bundle_out = bundle;

	return true;
}

typedef set<pair<int, int> > IntronSet;
typedef map<RefID, IntronSet> IntronTable;

struct AlignmentStats
{
    AlignmentStats()
    {
        _target_read_alignments = 0;
        _ref_read_alignments = 0;
        
        _fp_read_alignments = 0;
        _tp_read_alignments = 0;
        _fn_read_alignments = 0;
    }
    
    void register_false_read(const ReadHit& hit)
    {
        _target_read_alignments++;
        _fp_read_alignments++;
        
        
        // For now, we'll add these gaps to FP introns, but some of them
        // might actually be real, so we'll need to go through the list of TPs
        // when the user asks for the current FP rate.
        vector<pair<int, int> > introns; 
        hit.gaps(introns);
        
        pair<IntronTable::iterator, bool> ret = 
            _fp_introns.insert(make_pair(hit.ref_id(), IntronSet()));
        IntronSet& s = ret.first->second;
        
        copy(introns.begin(), introns.end(), inserter(s, s.end()));
    }
    
    void register_missed_read(const ReadHit& hit)
    {
        _ref_read_alignments++;
        _fn_read_alignments++;
        
        // For now, we'll add these gaps to FN introns, but some of them
        // might actually be real, so we'll need to go through the list of TPs
        // when the user asks for the current FN rate.
        vector<pair<int, int> > introns; 
        hit.gaps(introns);
        
        pair<IntronTable::iterator, bool> ret = 
            _fn_introns.insert(make_pair(hit.ref_id(), IntronSet()));
        IntronSet& s = ret.first->second;
        
        copy(introns.begin(), introns.end(), inserter(s, s.end()));  
    }
    
    void register_true_read(const ReadHit& hit)
    {
        vector<pair<int, int> > introns; 
        hit.gaps(introns);
        
        pair<IntronTable::iterator, bool> ret = 
            _tp_introns.insert(make_pair(hit.ref_id(), IntronSet()));
        IntronSet& s = ret.first->second;
        
        copy(introns.begin(), introns.end(), inserter(s, s.end())); 
        
        _target_read_alignments++;
        _ref_read_alignments++;
        _tp_read_alignments++;
    }
    
    int target_read_alignments() const { return _target_read_alignments; }
    int ref_read_alignments() const { return _ref_read_alignments; }
    int fp_read_alignments() const { return _fp_read_alignments; }
    int tp_read_alignments() const { return _tp_read_alignments; }
    int fn_read_alignments() const { return _fn_read_alignments; }
    
    int fp_introns() const 
    { 
        int num_fp_introns = 0;
        
        IntronTable::const_iterator fp_itr = _fp_introns.begin();
        while (fp_itr != _fp_introns.end())
        {
            const pair<RefID, IntronSet>& fps = *fp_itr;
            IntronTable::const_iterator ref_itr;
            ref_itr = _tp_introns.find(fps.first);
            if (ref_itr != _tp_introns.end())
            {
                IntronSet genuine_fps;
                set_difference(fps.second.begin(), 
                               fps.second.end(), 
                               ref_itr->second.begin(), 
                               ref_itr->second.end(),
                               inserter(genuine_fps, genuine_fps.begin()));
                num_fp_introns += genuine_fps.size();
            }
            ++fp_itr;
        }
        
        return num_fp_introns; 
    }
    
    int tp_introns() const 
    { 
        return (int)_tp_introns.size(); 
    }
    
    int fn_introns() const 
    { 
        int num_fn_introns = 0;
        
        IntronTable::const_iterator fn_itr = _fn_introns.begin();
        while (fn_itr != _fn_introns.end())
        {
            const pair<RefID, IntronSet>& fns = *fn_itr;
            IntronTable::const_iterator ref_itr;
            ref_itr = _tp_introns.find(fns.first);
            if (ref_itr != _tp_introns.end())
            {
                IntronSet genuine_fns;
                set_difference(fns.second.begin(), 
                               fns.second.end(), 
                               ref_itr->second.begin(), 
                               ref_itr->second.end(),
                               inserter(genuine_fns, genuine_fns.begin()));
                num_fn_introns += genuine_fns.size();
            }
            ++fn_itr;
        }
        
        return num_fn_introns; 
    }
    
private:
    int _target_read_alignments;
    int _ref_read_alignments;
    int _fp_read_alignments;
    int _tp_read_alignments;
    int _fn_read_alignments;
    
    int ref_introns;
    int target_introns;
    
    IntronTable _fp_introns;
    IntronTable _tp_introns;
    IntronTable _fn_introns;
};

void register_missed_read_alignment(const ReadHit& hit, 
                                    AlignmentStats& stats,
                                    FILE* fout)
{
    stats.register_missed_read(hit);
    
    if (fout)
    {
        //const string& buf = hit.orig_hit_str();
        fprintf(fout, "%s\n", hit.orig_hit_str().c_str());
    }
}

void register_missed_fragment_alignment(const MateHit& hit, 
                                        AlignmentStats& stats,
                                        FILE* fout)
{
    if (hit.left_alignment())
    {
        register_missed_read_alignment(*(hit.left_alignment()), stats, fout);
    }
    if (hit.right_alignment())
    {
        register_missed_read_alignment(*(hit.right_alignment()), stats, fout);
    }
}


void register_false_read_alignment(const ReadHit& hit, 
                                   AlignmentStats& stats,
                                   FILE* fout)
{
    stats.register_false_read(hit);
    
    if (fout)
    {
        //const string& buf = hit.orig_hit_str();
        fprintf(fout, "%s\n", hit.orig_hit_str().c_str());
    }
}


void register_false_fragment_alignment(const MateHit& hit, 
                                       AlignmentStats& stats,
                                       FILE* fout)
{
    if (hit.left_alignment())
    {
        register_false_read_alignment(*(hit.left_alignment()), stats, fout);
    }
    if (hit.right_alignment())
    {
        register_false_read_alignment(*(hit.right_alignment()), stats, fout);
    }
}

void register_true_read_alignment(const ReadHit& hit, 
                                  AlignmentStats& stats,
								  FILE* fout)
{

    stats.register_true_read(hit);
	
	if (fout)
    {
        //const string& buf = hit.orig_hit_str();
        fprintf(fout, "%s\n", hit.orig_hit_str().c_str());
    }
}

void register_true_fragment_alignment(const MateHit& hit, 
                                      AlignmentStats& stats,
									  FILE* fout)
{
    if (hit.left_alignment())
    {
        register_true_read_alignment(*(hit.left_alignment()), stats, fout);
    }
    if (hit.right_alignment())
    {
        register_true_read_alignment(*(hit.right_alignment()), stats, fout);
    }
}

void compare_bundles(const vector<MateHit>& ref_hits,
                     const vector<MateHit>& target_hits,
                     AlignmentStats& alignment_stats,
                     FILE* fmissed_hits,
                     FILE* ffalse_hits,
					 FILE* fcorrect_hits)
{
    vector<MateHit> true_positives;
    vector<MateHit> missed_ref_hits;
    vector<MateHit> false_target_hits;
    
    set_intersection(ref_hits.begin(), 
                     ref_hits.end(),
                     target_hits.begin(), 
                     target_hits.end(),
                     back_inserter(true_positives));
    
    set_difference(ref_hits.begin(), 
                   ref_hits.end(), 
                   target_hits.begin(), 
                   target_hits.end(),
                   back_inserter(missed_ref_hits));
    
    set_difference(target_hits.begin(), 
                   target_hits.end(), 
                   ref_hits.begin(), 
                   ref_hits.end(),
                   back_inserter(false_target_hits));
    
    
    foreach (const MateHit& m, false_target_hits)
    {
        register_false_fragment_alignment(m, alignment_stats, ffalse_hits);
    }
    
    foreach (const MateHit& m, missed_ref_hits)
    {
        register_missed_fragment_alignment(m, alignment_stats, fmissed_hits);
    }
    
    foreach (const MateHit& m, true_positives)
    {
        register_true_fragment_alignment(m, alignment_stats, fcorrect_hits);
    }
}

void print_alignment_stats(FILE* fout, const AlignmentStats& stats)
{
    double read_positives = stats.tp_read_alignments() + stats.fn_read_alignments();
    double read_guesses = stats.tp_read_alignments() + stats.fp_read_alignments();
    
    double read_recall = stats.tp_read_alignments() / read_positives;
    double read_precision = stats.tp_read_alignments() / read_guesses;
    
    double intron_positives = stats.tp_introns() + stats.fn_introns();
    double intron_guesses = stats.tp_introns() + stats.fp_introns();
    
    double intron_recall = stats.tp_introns() / intron_positives;
    double intron_precision = stats.tp_introns() / intron_guesses;
    
    fprintf(fout, "Total read alignment precision\t%lf\n", read_precision);
    fprintf(fout, "Total read alignment recall\t%lf\n", read_recall);
    
    fprintf(fout, "Intron precision\t%lf\n", intron_precision);
    fprintf(fout, "Intron recall\t%lf\n", intron_recall);
}

void get_ref_names(FILE* sam_in, set<string>& ref_names)
{
    char bwt_buf[2048];
	while(fgets(bwt_buf, 2048, sam_in))
    {
        // Are we still in the header region?
        if (bwt_buf[0] == '@')
            continue;
        
        const char* buf = bwt_buf;
        const char* _name = strsep((char**)&buf,"\t");
        if (!_name)
            continue;
        char name[2048];
        strncpy(name, _name, 2047); 
        
        const char* sam_flag_str = strsep((char**)&buf,"\t");
        if (!sam_flag_str)
           continue;
        
        const char* text_name = strsep((char**)&buf,"\t");
        if (!text_name)
            continue;
        ref_names.insert(text_name);
    }
    
}

void driver(FILE* ref_sam, 
            FILE* target_sam, 
            FILE* fmissed_hits, 
            FILE* ffalse_hits,
			FILE* fcorrect_hits)
{
    RefSequenceTable rt(true, false);
    
    set<string> ref_names;
    get_ref_names(ref_sam, ref_names);
    get_ref_names(target_sam, ref_names);
    foreach (const string& name, ref_names)
    {
        rt.get_id(name, NULL);
    }
    
    rewind(ref_sam);
    rewind(target_sam);
    
    ReadTable it;
    SAMHitFactory hs(it, rt);
    PositionBundleFactory ref_sam_factory(hs, ref_sam);
    PositionBundleFactory target_sam_factory(hs, target_sam);
    
    HitBundle curr_ref_bundle;
    HitBundle curr_target_bundle;
    
    bool valid_ref = ref_sam_factory.next_bundle(curr_ref_bundle);
    RefID ref_chr_order = rt.observation_order(curr_ref_bundle.ref_id());
    int ref_left = curr_ref_bundle.left();
    
    bool valid_target = target_sam_factory.next_bundle(curr_target_bundle);
    RefID targ_chr_order = rt.observation_order(curr_target_bundle.ref_id());
    int targ_left = curr_target_bundle.left();
   
    AlignmentStats alignment_stats;
    
    while (valid_target && valid_ref)
    {
        bool advance_target = false;
        bool advance_ref = false;
        
        if (ref_chr_order == targ_chr_order)
        {
            // compare the alignment coords (and maybe the CIGARs too)
            if (ref_left < targ_left)
            {
                // Then we missed this reference position => drop in recall.
                // need to advance the curr_ref_bundle
                foreach (const MateHit& hit, curr_ref_bundle.hits())
                {
                    register_missed_fragment_alignment(hit, 
                                                       alignment_stats,
                                                       fmissed_hits);
                }
                advance_ref = true;
                curr_ref_bundle = HitBundle();
            }
            else if (ref_left > targ_left)
            {
                // The hits in curr_target_bundle
                // are all false positives => drop in precision
                foreach (const MateHit& hit, curr_target_bundle.hits())
                {
                    register_false_fragment_alignment(hit, 
                                                      alignment_stats, 
                                                      ffalse_hits);
                }
                advance_target = true;
                curr_target_bundle = HitBundle();
            }
            else
            {
                //compare the reads in the two bundles
                compare_bundles(curr_ref_bundle.hits(),
                                curr_target_bundle.hits(),
                                alignment_stats,
                                fmissed_hits,
                                ffalse_hits,
								fcorrect_hits);
                advance_target = true;
                advance_ref = true;
                curr_target_bundle = HitBundle();
                curr_ref_bundle = HitBundle();
            }
        }
        else if (ref_chr_order < targ_chr_order)
        {
            // Then we missed this reference bundle => drop in recall.
            // need to advance the curr_ref_bundle
            foreach (const MateHit& hit, curr_ref_bundle.hits())
            {
                register_missed_fragment_alignment(hit, 
                                                   alignment_stats,
                                                   fmissed_hits);
            }
            curr_ref_bundle = HitBundle();
            advance_ref = true;
        }
        else 
        {
            // ref_chr_order > targ_chr_order. The hits in curr_target_bundle
            // are all false positives, and thus this case hurts precision
            foreach (const MateHit& hit, curr_target_bundle.hits())
            {
                register_false_fragment_alignment(hit, 
                                                  alignment_stats,
                                                  ffalse_hits);
            }
            curr_target_bundle = HitBundle();
            advance_target = true;
        }
        
        if (advance_target)
        {
            curr_target_bundle = HitBundle();
            
            valid_target = target_sam_factory.next_bundle(curr_target_bundle);
            targ_chr_order = rt.observation_order(curr_target_bundle.ref_id());
            targ_left = curr_target_bundle.left();
        }
        
        if (advance_ref)
        {
            curr_ref_bundle = HitBundle();
            valid_ref = ref_sam_factory.next_bundle(curr_ref_bundle);
            ref_chr_order = rt.observation_order(curr_ref_bundle.ref_id());
            ref_left = curr_ref_bundle.left();
        }
    } 
    
    while (valid_target)
    {
        foreach (const MateHit& hit, curr_target_bundle.hits())
        {
            register_false_fragment_alignment(hit, 
                                              alignment_stats,
                                              ffalse_hits);
        }
        curr_target_bundle = HitBundle();
        valid_target = ref_sam_factory.next_bundle(curr_target_bundle);
    }
    
    while (valid_ref)
    {
        foreach (const MateHit& hit, curr_ref_bundle.hits())
        {
            register_missed_fragment_alignment(hit, 
                                               alignment_stats,
                                               fmissed_hits);
        }
        curr_ref_bundle = HitBundle();
        valid_ref = ref_sam_factory.next_bundle(curr_ref_bundle);
    }
    
    print_alignment_stats(stderr, alignment_stats);
    
}

int main(int argc, char** argv)
{
	int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;
    
    assert (ref_sam_name != "" && target_sam_name != "");
    
    FILE* ref_sam = fopen(ref_sam_name.c_str(), "r");
	if (!ref_sam)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
				ref_sam_name.c_str());
		exit(1);
	}
    
    FILE* target_sam = fopen(target_sam_name.c_str(), "r");
	if (!target_sam)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for reading\n",
				target_sam_name.c_str());
		exit(1);
	}
    
    
    FILE* fmissed_hits = fopen("missed.sam", "w");
	if (!fmissed_hits)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
				"missed.sam");
		exit(1);
	}
    
    FILE* ffalse_hits = fopen("false.sam", "w");
	if (!ffalse_hits)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
				"false.sam");
		exit(1);
	}

	FILE* fcorrect_hits = fopen("correct.sam", "w");
	if (!fcorrect_hits)
	{
		fprintf(stderr, "Error: cannot open SAM file %s for writing\n",
				"false.sam");
		exit(1);
	}
	
    driver(ref_sam, target_sam, fmissed_hits, ffalse_hits, fcorrect_hits);
}

