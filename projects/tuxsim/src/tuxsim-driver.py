#! /usr/bin/env python
# encoding: utf-8
"""
tuxsim-driver.py

Created by Cole Trapnell on 2010-4-11.
Copyright (c) 2010 Cole Trapnell. All rights reserved.
"""

import sys
try:
    import psyco
    psyco.full()
except ImportError:
    pass

import sys
import getopt
import subprocess
import errno
import os
import tempfile
import warnings
import shutil
import copy
from datetime import datetime, date, time
import ConfigParser
import math
import numpy
import numpy.random
from numpy.random import *
import bisect
from operator import itemgetter


## {{{ http://code.activestate.com/recipes/577197/ (r9)
from bisect import bisect_left, bisect_right

class SortedCollection(object):
    '''Sequence sorted by a key function.

    SortedCollection() is much easier to work with than using bisect() directly.
    It supports key functions like those use in sorted(), min(), and max().
    The result of the key function call is saved so that keys can be searched
    efficiently.

    Instead of returning an insertion-point which can be hard to interpret, the
    five find-methods return a specific item in the sequence. They can scan for
    exact matches, the last item less-than-or-equal to a key, or the first item
    greater-than-or-equal to a key.

    Once found, an item's ordinal position can be located with the index() method.
    New items can be added with the insert() and insert_right() methods.
    Old items can be deleted with the remove() method.

    The usual sequence methods are provided to support indexing, slicing,
    length lookup, clearing, copying, forward and reverse iteration, contains
    checking, item counts, item removal, and a nice looking repr.

    Finding and indexing are O(log n) operations while iteration and insertion
    are O(n).  The initial sort is O(n log n).

    The key function is stored in the 'key' attibute for easy introspection or
    so that you can assign a new key function (triggering an automatic re-sort).

    In short, the class was designed to handle all of the common use cases for
    bisect but with a simpler API and support for key functions.

    >>> from pprint import pprint
    >>> from operator import itemgetter

    >>> s = SortedCollection(key=itemgetter(2))
    >>> for record in [
    ...         ('roger', 'young', 30),
    ...         ('angela', 'jones', 28),
    ...         ('bill', 'smith', 22),
    ...         ('david', 'thomas', 32)]:
    ...     s.insert(record)

    >>> pprint(list(s))         # show records sorted by age
    [('bill', 'smith', 22),
     ('angela', 'jones', 28),
     ('roger', 'young', 30),
     ('david', 'thomas', 32)]

    >>> s.find_le(29)           # find oldest person aged 29 or younger
    ('angela', 'jones', 28)
    >>> s.find_lt(28)           # find oldest person under 28
    ('bill', 'smith', 22)
    >>> s.find_gt(28)           # find youngest person over 28
    ('roger', 'young', 30)

    >>> r = s.find_ge(32)       # find youngest person aged 32 or older
    >>> s.index(r)              # get the index of their record
    3
    >>> s[3]                    # fetch the record at that index
    ('david', 'thomas', 32)

    >>> s.key = itemgetter(0)   # now sort by first name
    >>> pprint(list(s))
    [('angela', 'jones', 28),
     ('bill', 'smith', 22),
     ('david', 'thomas', 32),
     ('roger', 'young', 30)]

    '''

    def __init__(self, iterable=(), key=None):
        self._given_key = key
        key = (lambda x: x) if key is None else key
        decorated = sorted((key(item), item) for item in iterable)
        self._keys = [k for k, item in decorated]
        self._items = [item for k, item in decorated]
        self._key = key

    def _getkey(self):
        return self._key

    def _setkey(self, key):
        if key is not self._key:
            self.__init__(self._items, key=key)

    def _delkey(self):
        self._setkey(None)

    key = property(_getkey, _setkey, _delkey, 'key function')

    def clear(self):
        self.__init__([], self._key)

    def copy(self):
        return self.__class__(self, self._key)

    def __len__(self):
        return len(self._items)

    def __getitem__(self, i):
        return self._items[i]

    def __iter__(self):
        return iter(self._items)

    def __reversed__(self):
        return reversed(self._items)

    def __repr__(self):
        return '%s(%r, key=%s)' % (
            self.__class__.__name__,
            self._items,
            getattr(self._given_key, '__name__', repr(self._given_key))
        )

    def __reduce__(self):
        return self.__class__, (self._items, self._given_key)

    def __contains__(self, item):
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return item in self._items[i:j]

    def index(self, item):
        'Find the position of an item.  Raise ValueError if not found.'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].index(item) + i

    def count(self, item):
        'Return number of occurrences of item'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        j = bisect_right(self._keys, k)
        return self._items[i:j].count(item)

    def insert(self, item):
        'Insert a new item.  If equal keys are found, add to the left'
        k = self._key(item)
        i = bisect_left(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def insert_right(self, item):
        'Insert a new item.  If equal keys are found, add to the right'
        k = self._key(item)
        i = bisect_right(self._keys, k)
        self._keys.insert(i, k)
        self._items.insert(i, item)

    def remove(self, item):
        'Remove first occurence of item.  Raise ValueError if not found'
        i = self.index(item)
        del self._keys[i]
        del self._items[i]

    def find(self, k):
        'Return first item with a key == k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i != len(self) and self._keys[i] == k:
            return self._items[i]
        raise ValueError('No item found with key equal to: %r' % (k,))

    def find_le(self, k):
        'Return last item with a key <= k.  Raise ValueError if not found.'
        i = bisect_right(self._keys, k)
        if i:
            return self._items[i-1]
        raise ValueError('No item found with key at or below: %r' % (k,))

    def find_lt(self, k):
        'Return last item with a key < k.  Raise ValueError if not found.'
        i = bisect_left(self._keys, k)
        if i:
            return self._items[i-1]
        raise ValueError('No item found with key below: %r' % (k,))

    def find_ge(self, k):
        'Return first item with a key >= equal to k.  Raise ValueError if not found'
        i = bisect_left(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key at or above: %r' % (k,))

    def find_gt(self, k):
        'Return first item with a key > k.  Raise ValueError if not found'
        i = bisect_right(self._keys, k)
        if i != len(self):
            return self._items[i]
        raise ValueError('No item found with key above: %r' % (k,))


use_message = '''
 tuxsim-driver is a driver for benchmarking Cuffdiff
 
 Usage:
     tuxsim-driver [options] 
     
 Options:
     -o/--output-dir                <string>    [ default: ./ ]
     -p/--num-threads               <int>       [ default: 1  ]
     --keep-tmp
     --version
     --help
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

output_dir = "./"
logging_dir = output_dir + "logs/"
run_log = None
run_cmd = None

tmp_dir = output_dir + "tmp/"
bin_dir = sys.path[0] + "/"
overdispersion = 1.0
dispersion_filename = None
extern_simexpr_filename = None
#ok_str = "\t\t\t\t[OK]\n"
fail_str = "\t[FAILED]\n"

diff_mode = "uniform"
no_deseq = False
num_replicates = 3
min_fold_change = 1.0
max_fold_change = 2.0
num_perturbed_genes = 0
min_frags_for_perturbation = 10
read_length = 0
frag_length_mean = 0
frag_length_std_dev = 0
experiment_type = "depth-curve"
deseq_mode = "intersection-strict"
class TestParams:
        
    class SystemParams:
        def __init__(self,
                     threads,
                     keep_tmp):
            self.threads = threads
            self.keep_tmp = keep_tmp
            
        def parse_options(self, opts):
            for option, value in opts:
                if option in ("-p", "--num-threads"):
                    self.threads = int(value)
                if option in ("--keep-tmp"):
                    self.keep_tmp = True
        
        def check(self):
            pass
    class TestParams:
        def __init__(self,
                   skip_tuxsim,
                   skip_tophat,
                   skip_tophat_accuracy,
                   skip_asm,
                   skip_asm_perfect,
                   skip_asm_tophat,
                   skip_abund,
                   skip_abund_perfect,
                   skip_abund_tophat,
                   skip_replicate_diff):
            self.skip_tuxsim = skip_tuxsim
            self.skip_tophat = skip_tophat
            self.skip_tophat_accuracy = skip_tophat_accuracy
            self.skip_asm = skip_asm  
            self.skip_asm_perfect = skip_asm_perfect
            self.skip_asm_tophat = skip_asm_tophat
            self.skip_abund = skip_abund
            self.skip_abund_perfect = skip_abund_perfect
            self.skip_abund_tophat = skip_abund_tophat
            self.skip_replicate_diff = skip_replicate_diff
            
        def parse_options(self, opts):
            for option, value in opts:
                if option in ("--skip-tuxsim"):
                    self.skip_tuxsim = True
                if option in ("--skip-tophat"):
                    self.skip_tophat = True
                if option in ("--skip-tophat-accuracy"):
                    self.skip_tophat_accuracy = True
                if option in ("--skip-assemble"):
                    self.skip_asm = True
                if option in ("--skip-assemble-perfect-map"):
                    self.skip_asm_perfect = True
                if option in ("--skip-assemble-tophat-map"):
                    self.skip_asm_tophat = True
                if option in ("--skip-assemble"):
                    self.skip_abund = True
                if option in ("--skip-abundance-perfect-map"):
                     self.skip_abund_perfect = True
                if option in ("--skip-abundance-tophat-map"):
                    self.skip_abund_tophat = True
                if option in ("--skip-replicate-diff"):
                    self.skip_replicate_diff = True

        def check(self):
            pass
                   
    def __init__(self):        
        self.system_params = self.SystemParams(1,               # threads
                                               False)           # keep_tmp
        self.test_params = self.TestParams(False,   # skip_tuxsim
                                           False,   # skip_tophat
                                           False,   # skip_tophat_accuracy,
                                           False,   # skip_asm
                                           False,   # skip_asm_pefect
                                           False,   # skip_asm_tophat
                                           False,   # skip_abund
                                           False,   # skip_abund_perfect
                                           False,   # skip_abund_tophat
                                           False)   # skip_replicate_diff
    
    
    def check(self):
        self.system_params.check()
                   
    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvp:o:", 
                                        ["version",
                                         "help",  
                                         "output-dir=",
                                         "num-threads=",
                                         "keep-tmp",
                                         "overdispersion=",
                                         "dispersion-file=",
                                         "expression-file=",
                                         "num-replicates=",
                                         "min-fold-change=",
                                         "max-fold-change=",
                                         "diff-mode=",
                                         "num-genes=",
                                         "min-frags-for-perturb=",
                                         "skip-tuxsim",
                                         "skip-tophat",
                                         "skip-tophat-accuracy",
                                         "skip-assemble",
                                         "skip-assemble-perfect-map",
                                         "skip-assemble-tophat-map",
                                         "skip-abundance",
                                         "skip-abundance-perfect-map",
                                         "skip-abundance-tophat-map",
                                         "skip-replicate-diff",
                                         "experiment-type=",
                                         "deseq-mode="])
        except getopt.error, msg:
            raise Usage(msg)
            
        self.system_params.parse_options(opts)
        self.test_params.parse_options(opts)
       
        # option processing
        for option, value in opts:
            if option in ("-v", "--version"):
                print "tuxedo_test v%s" % (get_version())
                exit(0)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option == "--overdispersion":
                global overdispersion
                overdispersion = float(value)
            if option == "--dispersion-file":
                global dispersion_filename
                dispersion_filename = value
            if option == "--expression-file":
                global extern_simexpr_filename
                extern_simexpr_filename = value
            if option == "--diff-mode":
                global diff_mode
                diff_mode = value
            if option == "--num-replicates":
                global num_replicates
                num_replicates = int(value)
            if option == "--max-fold-change":
                global max_fold_change
                max_fold_change = float(value)
            if option == "--min-fold-change":
                global min_fold_change
                min_fold_change = float(value)
            # if option == "--no-deseq":
            #     global no_deseq
            #     no_deseq = True
            if option == "--experiment-type":
                global experiment_type
                experiment_type = value
            if option == "--deseq-mode":
                global deseq_mode
                deseq_mode = value
            if option == "--num-genes":
                global num_perturbed_genes
                num_perturbed_genes = int(value)
                #print >> sys.stderr, "NUM GENES = %d" % num_perturbed_genes
            if option == "--min-frags-for-perturb":
                global min_frags_for_perturbation
                min_frags_for_perturbation = int(value)
            if option in ("-o", "--output-dir"):
                global output_dir
                global logging_dir
                global tmp_dir
                output_dir = value + "/"
                logging_dir = output_dir + "logs/"
                tmp_dir = output_dir + "tmp/"
            
        return args
    

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

def prepare_output_dir():
    
    print >> sys.stderr, "[%s] Preparing output location %s" % (right_now(), output_dir)
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)
        
    if os.path.exists(tmp_dir):
        pass
    else:        
        os.mkdir(tmp_dir)
        
def formatTD(td):
    hours = td.seconds // 3600
    minutes = (td.seconds % 3600) // 60
    seconds = td.seconds % 60
    return '%02d:%02d:%02d' % (hours, minutes, seconds) 

def tmp_name():
    tmp_root = output_dir + "tmp/"
    if os.path.exists(tmp_root):
        pass
    else:        
        os.mkdir(tmp_root)
    return tmp_root + os.tmpnam().split('/')[-1] 

def get_version():
   return "0.0.1"

def tuxsim(config_file, params, external_abundance_filename=None, out_prefix=None):
    print >> sys.stderr, "[%s] Generating RNA-Seq reads" % (right_now())

    cmd = ["tuxsim"]
    
    if external_abundance_filename != None:
        cmd.append("--expression")
        cmd.append(external_abundance_filename)
    if out_prefix != None:
        cmd.append("--output.prefix")
        cmd.append(out_prefix)
        
    cmd.append(config_file)

    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute tuxsim"
            exit(1)
    # tuxsim not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: tuxsim not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
        
def tophat(params):
    print >> sys.stderr, "[%s] Aligning RNA-Seq reads" % (right_now())

    cmd = ["tophat"]
    
    # Run TopHat with more than one thread?
    cmd.extend(["-p", str(params.system_params.threads)])
    cmd.extend(["-F 0.0"])

    # TODO: pick up library fragment mean from sim.cfg instead of using 
    # hardcoded value.
    cmd.extend(["-r","50"])
    
    # Add required Bowtie index and read args to TopHat command
    cmd.extend(["index/ref", "sim_1.fq", "sim_2.fq"])
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        # Tophat reported an error
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute tophat"
            exit(1)

        # Add required Bowtie index and read args to TopHat command
        convert_cmd = ["samtools", "view", "tophat_out/accepted_hits.bam"]
        try:       
            print >> run_log, " ".join(convert_cmd)
            ret = subprocess.call(convert_cmd, stdout=open("tophat_out/accepted_hits.sam", "w"))
                                  
            if ret != 0:
                print >> sys.stderr, fail_str, "Error: could not execute tophat"
                exit(1)
        except OSError, o:
            if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
                print >> sys.stderr, fail_str, "Error: samtools not found on this system.  Did you forget to include it in your PATH?"
            exit(1)

        
    # tophat not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: tophat not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
        
def cufflinks(params, out_dir, sam_file, gtf_file=None):
    if gtf_file != None:
        print >> sys.stderr, "[%s] Quantitating transcripts" % (right_now())
    else:
        print >> sys.stderr, "[%s] Assembling transcripts" % (right_now())

    cmd = ["cufflinks"]
    
    if out_dir != None and out_dir != "":
        cmd.extend(["-o", out_dir])
    
    if gtf_file != None:
        cmd.extend(["-G", gtf_file])
    
    # Run Cufflinks with more than one thread?
    cmd.extend(["-p", str(params.system_params.threads)])
        
    cmd.append(sam_file)
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cufflinks"
            exit(1)
    # cufflinks not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cufflinks not found on this system.  Did you forget to include it in your PATH?"
        exit(1) 
        
def cuffdiff(params, 
             out_dir, 
             gtf_file, 
             sam_file1, 
             sam_file2, 
             labels=None, 
             min_frags=None, 
             read_skip_fraction=None,
             no_read_pairing=False,
             read_trim_length=None,
             library_type="fr-unstranded",
             fdr=0.01):

    print >> sys.stderr, "[%s] Testing for differences with cuffdiff" % (right_now())

    cmd = ["cuffdiff"]

    if out_dir != None and out_dir != "":
        cmd.extend(["-o", out_dir])

    # Run Cuffdiff with more than one thread?
    cmd.extend(["-p", str(params.system_params.threads)])
    
    if min_frags != None:
        cmd.extend(["-c", str(min_frags)])
    if labels != None:
        cmd.extend(["-L", ",".join(labels)])
    cmd.extend(["--emit-count-tables"])  
    #cmd.extend(["--emit-count-tables"])  
            
    if read_skip_fraction != None:
        cmd.extend(["--read-skip-fraction", str(read_skip_fraction)])
    
    if read_trim_length != None:
        cmd.extend(["--trim-read-length", str(read_trim_length)]) 
        
    if no_read_pairing == True:
        cmd.extend(["--no-read-pairs"])
    
    global frag_length_mean    
    cmd.extend(["-m", str(frag_length_mean)])
        
    global frag_length_std_dev 
    cmd.extend(["-s", str(frag_length_std_dev)])
    
    cmd.extend(["--FDR", str(fdr)])
    cmd.extend(["--library-type", library_type])
            
    cmd.append(gtf_file)
    cmd.append(sam_file1)
    cmd.append(sam_file2)

    
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffdiff"
            exit(1)
    # cuffdiff not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffdiff not found on this system.  Did you forget to include it in your PATH?"
        exit(1)


def samcompare(params, out_dir, ref_sam, target_sam, report_sam_files=True):
    
    print >> sys.stderr, "[%s] Assessing fragment alignment accuracy" % (right_now())
        
    cmd = ["samcompare"]
    if out_dir != None and out_dir != "":
        cmd.extend(["--output-dir", out_dir])
    if report_sam_files == True:
        cmd.extend(["--report-sam"])
    
    cmd.append(ref_sam)
    cmd.append(target_sam)
    
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute samcompare"
            exit(1)
    # samcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: samcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)       
        
def cuffcompare(params, prefix, ref_gtf, cuff_gtf):
    
    print >> sys.stderr, "[%s] Comparing reference %s to assembly %s" % (right_now(), ref_gtf, cuff_gtf)
    cmd = ["cuffcompare"]

    cmd.extend(["-o", prefix])
    cmd.extend(["-r", ref_gtf])
    cmd.extend([cuff_gtf])
    
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute cuffcompare"
            exit(1)
    # cuffcompare not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: cuffcompare not found on this system.  Did you forget to include it in your PATH?"
        exit(1)


# Perturbes the given expr file by permuting the rhos, leaving the 
# total expression for each gene the same 
def create_rotate_rhos_condition(original_sample_filename, 
                                      cond_simexpr,
                                      cond_perturbed,
                                      num_perturbed_genes,
                                      total_frags):

    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)

    original_sample = open(original_sample_filename, "r")

    all_genes = {}

    total_fpkm = 0.0
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        total_fpkm += fpkm
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append([trans_id, rho, eff_len])

    expressed_genes = set([])
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2

        if fpkm > 0 and ((rho * total_frags) / (eff_len / 1000.0)) > min_frags_for_perturbation and len (all_isos_for_gene) >= min_isos_for_perturb:
            trans_list = expressed_genes.add(gene_id)

    #print len(expressed_genes)
    #expressed_genes = list(expressed_genes)

    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)

    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return

    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")

    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    perturbed_gene_ids = set([])

    new_rho_for_trans_ids = {}
    max_try = 0

    expressed_genes = list(expressed_genes)
    #print expressed_genes
    while len(perturbed_gene_ids) < num_perturbed_genes and max_try < 10000 * len(expressed_genes):
        #trans_ids = []
        p_i = randint(0, len(expressed_genes)-1)
        max_try += 1

        # if (max_try % 100 == 0):
        #             print max_try
        if expressed_genes[p_i] not in perturbed_gene_ids:
            p = expressed_genes[p_i]
        else:
            continue
        id_list = all_genes.get(p)

        id_list.sort(lambda x,y: cmp(y[2],x[2]))

        counts = [r * total_frags for [trans_id, r, length] in id_list]
        rhos = [r for [trans_id, r,  length] in id_list]
        lengths = [length for [trans_id, r, length] in id_list]

        new_rhos = list(rhos)
        shuffle(new_rhos)
        
        #print "***********"
        #print rhos, sum(rhos)
        #print new_rhos, sum(new_rhos)

        # print "***********"
        # print p
        # print rhos, sum(rhos)
        # print new_rhos, sum(new_rhos)
        for i in range(0, len(counts)):
            id_list[i][1] = new_rhos[i]
        perturbed_gene_ids.add(p)
        all_genes[p] = id_list

    print >> sys.stderr, "After %d attempts, only perturbed %d genes (of %d expressed)" % (max_try, len(perturbed_gene_ids), len(expressed_genes))
    transcript_num = 0
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1

        if len(cols) < 7:
            continue
        gene_id = cols[0]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])

        if gene_id in perturbed_gene_ids:
            #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * max_fold_change
            #rho *= max_fold_change
            id_list = all_genes[gene_id]
            for (trans_id, r, eff_len) in id_list:
                if trans_id == cols[1]:
                    rho = r
        
        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs

    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 

    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[0] in perturbed_gene_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()
    
# Perturbes the given expr file by reflecting the rhos, leaving the 
# total expression for each gene the same 
def create_reflect_isoforms_condition(original_sample_filename, 
                                      cond_simexpr,
                                      cond_perturbed,
                                      num_perturbed_genes,
                                      total_frags):

    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)

    original_sample = open(original_sample_filename, "r")

    all_genes = {}

    total_fpkm = 0.0
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        total_fpkm += fpkm
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append([trans_id, rho, eff_len])

    expressed_genes = set([])
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2

        if fpkm > 0 and ((rho * total_frags) / (eff_len / 1000.0)) > min_frags_for_perturbation and len (all_isos_for_gene) >= min_isos_for_perturb:
            trans_list = expressed_genes.add(gene_id)

    #print len(expressed_genes)
    #expressed_genes = list(expressed_genes)

    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)

    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return

    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")

    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    perturbed_gene_ids = set([])

    new_rho_for_trans_ids = {}
    max_try = 0

    expressed_genes = list(expressed_genes)
    #print expressed_genes
    while len(perturbed_gene_ids) < num_perturbed_genes and max_try < 10000 * len(expressed_genes):
        #trans_ids = []
        p_i = randint(0, len(expressed_genes)-1)
        max_try += 1

        # if (max_try % 100 == 0):
        #             print max_try
        if expressed_genes[p_i] not in perturbed_gene_ids:
            p = expressed_genes[p_i]
        else:
            continue
        id_list = all_genes.get(p)

        id_list.sort(lambda x,y: cmp(x[1],y[1]))
        
        least_abundant_idx = 0
        least_abundant = 0.0
        most_abundant = id_list[-1][1]
        while least_abundant_idx != len(id_list):
            if id_list[least_abundant_idx][1] != 0:
                least_abundant = id_list[least_abundant_idx][1]
                break
            least_abundant_idx += 1
        
        if least_abundant == 0 or least_abundant == most_abundant:
            continue
        
        # Only select this gene if the major and minor isoforms are actually
        # different.
        try:
            if abs(math.log(most_abundant/least_abundant, 2)) < 1.0:
                continue
        except ValueError, e:
            print most_abundant, least_abundant
        counts = [r * total_frags for [trans_id, r, length] in id_list]
        rhos = [r for [trans_id, r,  length] in id_list]
        lengths = [length for [trans_id, r, length] in id_list]
        
        new_rhos = list(rhos)
        #shuffle(new_rhos)
        new_rhos.reverse()

        #print "***********"
        #print rhos, sum(rhos)
        #print new_rhos, sum(new_rhos)

        # print "***********"
        # print p
        # print rhos, sum(rhos)
        # print new_rhos, sum(new_rhos)
        for i in range(0, len(counts)):
            id_list[i][1] = new_rhos[i]
        perturbed_gene_ids.add(p)
        all_genes[p] = id_list

    print >> sys.stderr, "After %d attempts, only perturbed %d genes (of %d expressed)" % (max_try, len(perturbed_gene_ids), len(expressed_genes))
    transcript_num = 0
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1

        if len(cols) < 7:
            continue
        gene_id = cols[0]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])

        if gene_id in perturbed_gene_ids:
            #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * max_fold_change
            #rho *= max_fold_change
            id_list = all_genes[gene_id]
            for (trans_id, r, eff_len) in id_list:
                if trans_id == cols[1]:
                    rho = r

        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs

    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 

    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[0] in perturbed_gene_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()

# Perturbes the given expr file by shifting all of selected genes' splice
# variancts onto the shortest splice variance of each
def create_shift_all_to_one_condition(original_sample_filename, 
                                      cond_simexpr,
                                      cond_perturbed,
                                      max_fold_change,
                                      num_perturbed_genes,
                                      total_frags):

    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)

    original_sample = open(original_sample_filename, "r")

    all_genes = {}

    total_fpkm = 0.0
    total_count = 0.0
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        count = fpkm * (eff_len/1000.0) * (total_frags/1000000.0)
        total_count += count
        total_fpkm += fpkm
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append([trans_id, rho, fpkm, count, eff_len])

    expressed_genes = set([])
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2
        
        if (fpkm > 0) and ((rho * total_frags) / (eff_len / 1000.0) > min_frags_for_perturbation) and len (all_isos_for_gene) >= min_isos_for_perturb:
            trans_list = expressed_genes.add(gene_id)

    #print len(expressed_genes)
    #expressed_genes = list(expressed_genes)

    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)

    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return
    
    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")

    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    perturbed_gene_ids = set([])

    new_rho_for_trans_ids = {}
    max_try = 0

    expressed_genes = list(expressed_genes)
    #print expressed_genes
    while len(perturbed_gene_ids) < num_perturbed_genes and max_try < 100 * len(expressed_genes):
        #trans_ids = []
        p_i = randint(0, len(expressed_genes)-1)
        max_try += 1

        # if (max_try % 100 == 0):
        #             print max_try
        if expressed_genes[p_i] not in perturbed_gene_ids:
            p = expressed_genes[p_i]
        else:
            continue
        id_list = all_genes.get(p)

        #id_list.sort(lambda x,y: cmp(y[2],x[2]))

        counts = [count for [trans_id, r, fpkm, count, length] in id_list]
        rhos = [r for [trans_id, r, fpkm, count, length] in id_list]
        lengths = [length for [trans_id, r, fpkm, count, length] in id_list]
        
        new_rhos = []
        new_iso = min(lengths)
        shifted = False
        print "###############"
        for i in range(0, len(rhos)):
            if shifted == False and lengths[i] == new_iso:
                new_fpkm = (sum(counts) / (new_iso/1000.0)) / (total_frags / 1000000.0)
                nr = new_fpkm / total_fpkm
                # print "transcript id = ", id_list[i][0]
                # print "old FPKM = ", id_list[i][2]
                # print "new FPKM = ", new_fpkm
                #print "old counts = %s (or %s)" % (sum(counts), sum(rhos) * total_count)
                #print "new counts = ", nr * total_count
                # print "length = ", lengths[i]
                # print "new rho = ", nr
                # print "old rho = ", rhos[i]
                # print "fold change in rho = ", nr/rhos[i]
                new_rhos.append(nr)
                shifted = True
            else:
                new_rhos.append(0)
        # 
        #print "***********"
        # print counts, sum(counts)
        # print rhos, sum(rhos)
        # print new_rhos, sum(new_rhos)
        # print lengths
        # print sum(new_rhos)/sum(rhos)
        if sum(rhos) > 0 and sum(new_rhos)/sum(rhos) > max_fold_change: 
            print "***********"
            print p
            print counts, sum(counts)
            print rhos, sum(rhos)
            print new_rhos, sum(new_rhos)
            print sum(new_rhos)/sum(rhos)
            for i in range(0, len(counts)):
                id_list[i][1] = new_rhos[i]
            perturbed_gene_ids.add(p)
            all_genes[p] = id_list
        else:
            #print "Insufficient fold change (%f, %f/%f), skipping" % (sum(new_rhos)/sum(rhos), sum(counts), new_iso)
            continue

    print >> sys.stderr, "After %d attempts, only perturbed %d genes (of %d expressed)" % (max_try, len(perturbed_gene_ids), len(expressed_genes))
    transcript_num = 0
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1

        if len(cols) < 7:
            continue
        gene_id = cols[0]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])

        if gene_id in perturbed_gene_ids:
            #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * max_fold_change
            #rho *= max_fold_change
            id_list = all_genes[gene_id]
            for (trans_id, r, fpkm, count, eff_len) in id_list:
                if trans_id == cols[1]:
                    rho = r

        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs

    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 

    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[0] in perturbed_gene_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()

# Creates a replicate, choosing a random set of genes and increasing or decreasing
# their expression by a given multiplier
def create_isoform_switch_condition(original_sample_filename, 
                                    cond_simexpr,
                                    cond_perturbed,
                                    max_fold_change,
                                    num_perturbed_genes,
                                    total_frags):

    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)

    original_sample = open(original_sample_filename, "r")

    all_genes = {}

    total_fpkm = 0.0
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        total_fpkm += fpkm
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append([trans_id, rho, eff_len])

    expressed_genes = set([])
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2]) 
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2
        if (fpkm > 0) and ((rho * total_frags) / (eff_len / 1000.0) > min_frags_for_perturbation) and len (all_isos_for_gene) > min_isos_for_perturb:
            trans_list = expressed_genes.add(gene_id)

    #print len(expressed_genes)
    #expressed_genes = list(expressed_genes)

    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)

    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return

    #perturbed_genes = [0 for i in range(0, num_expressed_genes)]
    # perturbed = numpy.random.random_integers(0, num_expressed_genes - 1, num_perturbed_genes)
    
    # expressed_gene_ids = list(expressed_genes)
    # for p in perturbed:
    #     perturbed_gene_ids.add(expressed_gene_ids[p])

    #print >> sys.stderr, "Perturbing %d genes by %f fold" % (len(perturbed_gene_ids), max_fold_change)
    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")

    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    perturbed_gene_ids = set([])

    new_rho_for_trans_ids = {}
    max_try = 0
    
    expressed_genes = list(expressed_genes)
    while len(perturbed_gene_ids) < num_perturbed_genes and max_try < 100 * len(expressed_genes):
        #trans_ids = []
        p_i = randint(0, len(expressed_genes)-1)
        max_try += 1
        
        #if (i % 100 == 0):
        #print max_try
        if expressed_genes[p_i] not in perturbed_gene_ids:
            p = expressed_genes[p_i]
        else:
            continue
        id_list = all_genes.get(p)
        
        
        id_list.sort(lambda x,y: cmp(y[2],x[2]))
        
        counts = [r * total_frags for [trans_id, r, length] in id_list]
        rhos = [r for [trans_id, r,  length] in id_list]
        lengths = [length for [trans_id, r, length] in id_list]
        #print rhos
        
        #shuffle(counts)
        total_gene_count = sum(counts)
        
        #new_rhos = [(counts[i]/(lengths[i]/1000.0)/(total_frags/1000000.0))/total_fpkm for i in range(0, len(counts))]
        
        #new_rhos = list(rhos)
        new_rhos = [2*x for x in rhos]
        shuffle(new_rhos)
        
        #print "***********"
        #print rhos, sum(rhos)
        #print new_rhos, sum(new_rhos)
        
        if sum(rhos) > 0 and sum(new_rhos)/sum(rhos) > max_fold_change: 
            # print "***********"
            # print p
            # print rhos, sum(rhos)
            # print new_rhos, sum(new_rhos)
            for i in range(0, len(counts)):
                id_list[i][1] = new_rhos[i]
            perturbed_gene_ids.add(p)
            all_genes[p] = id_list
        else:
            #print "Insufficient fold change, skipping"
            continue
    
    print >> sys.stderr, "After %d attempts, only perturbed %d genes (of %d expressed)" % (max_try, len(perturbed_gene_ids), len(expressed_genes))
    transcript_num = 0
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1

        if len(cols) < 7:
            continue
        gene_id = cols[0]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])

        if gene_id in perturbed_gene_ids:
            #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * max_fold_change
            #rho *= max_fold_change
            id_list = all_genes[gene_id]
            for (trans_id, r, eff_len) in id_list:
                if trans_id == cols[1]:
                    rho = r
         
        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs

    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 

    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[0] in perturbed_gene_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()

# Creates a replicate, choosing a random set of genes and increasing or decreasing
# their expression by a given multiplier
def create_perturbed_condition(original_sample_filename, 
                               cond_simexpr,
                               cond_perturbed,
                               max_fold_change, 
                               num_perturbed_genes,
                               total_frags):
    
    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)
    
    original_sample = open(original_sample_filename, "r")
    
    all_genes = {}
    
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append(trans_id)
    
    expressed_genes = {}
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            print >> sys.stderr, "Warning: malformed line during load for perturbation, skipping"
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2
        if fpkm > 0 and rho * total_frags / (eff_len / 1000.0) > min_frags_for_perturbation and len (all_isos_for_gene) > min_isos_for_perturb:
            trans_list = expressed_genes.setdefault(gene_id, [])
            trans_list.append(trans_id)
    
    #expressed_genes = list(expressed_genes)
    
    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)
    
    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return
    else:
        print >> sys.stderr, "Expressed set has %d candidates for perturbation" % num_expressed_genes
    
    #perturbed_genes = [0 for i in range(0, num_expressed_genes)]
    if num_expressed_genes > num_perturbed_genes and num_expressed_genes > 0:
        perturbed = numpy.random.random_integers(0, num_expressed_genes - 1, num_perturbed_genes)
    else:
        perturbed = []
    perturbed_gene_ids = set([])
    expressed_gene_ids = list(expressed_genes.keys())
    for p in perturbed:
        perturbed_gene_ids.add(expressed_gene_ids[p])
    
    print >> sys.stderr, "Perturbing %d genes by %f fold" % (len(perturbed_gene_ids), max_fold_change)
    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")
    
    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    total_fpkm = 0.0
    
    transcript_num = 0
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1
        
        if len(cols) < 7:
            print >> sys.stderr, "Warning: malformed line during perturbation, skipping"
            continue
        gene_id = cols[0]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])
                
        if gene_id in perturbed_gene_ids:
            #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * max_fold_change
            rho *= max_fold_change
        if rho < 0:
            print >> sys.stderr, "Warning: rho = %lg can't be negative, setting to 0.0" % rho
            rho = 0

        
        total_fpkm += float(cols[2]) 
        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs
    
    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 
    
    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[0] in perturbed_gene_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()

# Creates a replicate, choosing a random set of genes and increasing or decreasing
# their expression by a given multiplier
def create_interval_fold_change_condition(original_sample_filename, 
                                          cond_simexpr,
                                          cond_perturbed,
                                          min_fold_change,
                                          max_fold_change, 
                                          num_perturbed_genes,
                                          total_frags,
                                          iso_perturb_mode="single"):

    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)

    original_sample = open(original_sample_filename, "r")

    all_genes = {}

    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            print >> sys.stderr, "Error: bad in line interval fold change perturbation (1), skipping"
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append([trans_id, rho, eff_len])

    expressed_genes = {}
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            print >> sys.stderr, "Error: bad in line interval fold change perturbation (2), skipping"
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2
        #print >> sys.stderr, (rho * total_frags) / (eff_len / 1000.0)
        if fpkm > 0 and ((rho * total_frags) / (eff_len / 1000.0)) > min_frags_for_perturbation and len (all_isos_for_gene) >= min_isos_for_perturb:
            trans_list = expressed_genes.setdefault(gene_id, [])
            trans_list.append([trans_id, rho, eff_len])

    #expressed_genes = list(expressed_genes)

    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)

    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return
    else:
        print >> sys.stderr, "Expressed set has %d candidates for perturbation" % num_expressed_genes

    #perturbed_genes = [0 for i in range(0, num_expressed_genes)]
    if num_expressed_genes > num_perturbed_genes and num_expressed_genes > 0:
        perturbed = numpy.random.random_integers(0, num_expressed_genes - 1, num_perturbed_genes)
    else:
        perturbed = []
    perturbed_gene_ids = set([])
    expressed_gene_ids = list(expressed_genes.keys())
    for p in perturbed:
        perturbed_gene_ids.add(expressed_gene_ids[p])

    perturbed_trans_ids = set([])

    print >> sys.stderr, "Perturbing %d genes" % (len(perturbed_gene_ids))
    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")

    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    total_fpkm = 0.0

    transcript_num = 0
    
    # This shuffling procedure will let us pick a random isoform and upregulate it by 
    # a random amount.  That way, we aren't biasing ourselves only to the most highly expressed
    # isoforms
    for gene_id, trans_list in expressed_genes.iteritems():
        numpy.random.shuffle(trans_list)
    
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1

        if len(cols) < 7:
            print >> sys.stderr, "Error: bad in line interval fold change perturbation (3), skipping"
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])
        fold_change = numpy.random.uniform(min_fold_change, max_fold_change, 1)[0]
        direction = [-1,1]
        numpy.random.shuffle(direction)
        #print >> sys.stderr, direction
        fold_change *= direction[0]
        fold_change = 2.0**fold_change
        # Flip the top isoform up by fold_change
        if gene_id in perturbed_gene_ids:
            trans_list = expressed_genes.get(gene_id)
            #trans_list.sort(lambda x,y: cmp(y[1],x[1]))
            #print trans_list[0][0], trans_id
            if iso_perturb_mode == "single":
                # just perturb a random isoform (the list was shuffled above)
                if trans_id == trans_list[0][0]:
                    #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * fold_change
                    rho *= fold_change
                    perturbed_trans_ids.add(trans_id)
            else: 
                # perturb all the isoforms
                rho *= fold_change
                perturbed_trans_ids.add(trans_id)
                
        if rho < 0:
            print >> sys.stderr, "Warning: rho = %lg can't be negative, setting to 0.0" % rho
            rho = 0

        total_fpkm += float(cols[2]) 
        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs

    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 

    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[1] in perturbed_trans_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    print >> sys.stderr, "Wrote perturbed lines to %s" % cond_perturbed
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()
        
        
# Creates a replicate, choosing a random set of genes and increasing or decreasing
# their expression by a given multiplier
def create_single_iso_up_condition(original_sample_filename, 
                                   cond_simexpr,
                                   cond_perturbed,
                                   max_fold_change, 
                                   num_perturbed_genes,
                                   total_frags):

    print >> sys.stderr, "[%s] Creating a perturbation of %s" % (right_now(), original_sample_filename)

    original_sample = open(original_sample_filename, "r")

    all_genes = {}

    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        trans_list = all_genes.setdefault(gene_id, [])
        trans_list.append([trans_id, rho, eff_len])

    expressed_genes = {}
    original_sample.seek(0)
    original_sample.readline()
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        rho = float(cols[3])
        fpkm = float(cols[2])
        eff_len = float(cols[6])
        all_isos_for_gene = all_genes.get(gene_id)
        min_isos_for_perturb = 2
        #print >> sys.stderr, (rho * total_frags) / (eff_len / 1000.0)
        if fpkm > 0 and ((rho * total_frags) / (eff_len / 1000.0)) > min_frags_for_perturbation and len (all_isos_for_gene) >= min_isos_for_perturb:
            trans_list = expressed_genes.setdefault(gene_id, [])
            trans_list.append([trans_id, rho, eff_len])

    #expressed_genes = list(expressed_genes)

    num_expressed_genes = len(expressed_genes) - 1
    original_sample.seek(0)

    if num_expressed_genes < 0:
        print >> sys.stderr, "Error: no genes expressed in original sample!"
        return
    else:
        print >> sys.stderr, "Expressed set has %d candidates for perturbation" % num_expressed_genes

    #perturbed_genes = [0 for i in range(0, num_expressed_genes)]
    if num_expressed_genes > num_perturbed_genes and num_expressed_genes > 0:
        perturbed = numpy.random.random_integers(0, num_expressed_genes - 1, num_perturbed_genes)
    else:
        perturbed = []
    perturbed_gene_ids = set([])
    expressed_gene_ids = list(expressed_genes.keys())
    for p in perturbed:
        perturbed_gene_ids.add(expressed_gene_ids[p])
    
    perturbed_trans_ids = set([])
    
    print >> sys.stderr, "Perturbing %d genes by %f fold" % (len(perturbed_gene_ids), max_fold_change)
    rep = open(cond_simexpr, "w")
    perturbed_out = open(cond_perturbed, "w")

    header = original_sample.readline()

    total_rep_density = 0.0
    new_expr_recs = []
    total_fpkm = 0.0
    
    transcript_num = 0
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        #transcript_num += 1

        if len(cols) < 7:
            continue
        gene_id = cols[0]
        trans_id = cols[1]
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        rho = float(cols[3])

        # Flip the top isoform up by max_fold_change
        if gene_id in perturbed_gene_ids:
            trans_list = expressed_genes.get(gene_id)
            trans_list.sort(lambda x,y: cmp(y[1],x[1]))
            #print trans_list[0][0], trans_id
            if trans_id == trans_list[0][0]:
                #print >> sys.stderr, cols[0], frag_cov * eff_len, frag_cov * eff_len * max_fold_change
                rho *= max_fold_change
                perturbed_trans_ids.add(trans_id)

        if rho < 0:
            print >> sys.stderr, "Warning: rho = %lg can't be negative, setting to 0.0" % rho
            rho = 0
            
        total_fpkm += float(cols[2]) 
        total_rep_density += rho
        new_expr_recs.append([rho, cols])

    print >> rep, header,
    print >> perturbed_out, header,
    #print >> sys.stderr, new_expr_recs

    print >> sys.stderr, "Total perturbed rho =", total_rep_density 
    print >> sys.stderr, "Total perturbed fragments =", total_frags 

    rec_num = 0
    for rec in new_expr_recs:
        new_rho = rec[0] / total_rep_density
        cols = rec[1]
        eff_len = float(cols[6])
        new_fpkm = new_rho * total_fpkm
        if (eff_len > 0 and total_frags > 0) == False:
            new_fpkm = 0
        cols = rec[1]
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        if cols[1] in perturbed_trans_ids:
            print >> perturbed_out, "%s\t%s\t%s\t%g\t%s\t%s\t%s\t%s" % (cols[0], cols[1], new_fpkm, new_rho, cols[4], rec[0], cols[6], cols[7]) 
        rec_num += 1
    perturbed_out.flush()
    perturbed_out.close()
    rep.flush()
    rep.close()
        
def compress_gtf(params, ref_gtf, compressed_gtf, mode):

    print >> sys.stderr, "[%s] Compressing GTF %s into %s" % (right_now(), ref_gtf, compressed_gtf)
    cmd = ["compress_gtf"]
    
    cmd.extend([mode])
    cmd.extend([ref_gtf])
    cmd.extend([compressed_gtf])
    
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute compress_gtf"
            exit(1)
    # compress_gtf not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: compress_gtf not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
    
def load_empirical_dispersion(count_filename, allow_undispersion=False):
    count_f = open(count_filename)
    header = count_f.readline()

    s = SortedCollection(key=itemgetter(1))

    s.insert((0,0))
    
    total_frags = 0
    for line in count_f:
        line = line.strip()
        cols = line.split("\t")
        count = float(cols[0])
        var = float(cols[1])
        fitted_var = float(cols[2])
        total_frags += count
        if allow_undispersion == False:
            fitted_var = max(fitted_var, count)
        try:
            s.find(count)
            #print >> sys.stderr, "Table already contains value for %f" % count
        except ValueError, e:
            s.insert((count,fitted_var))

    return (s, total_frags)

def get_empirical_var(disp_table, count):
    try:
        lower = disp_table.find_le(count)
    except ValueError, e:
        print >> sys.stderr,  "WARNING: out of dynamic count range (too low)"
        exit(1)
    try:
        upper = disp_table.find_ge(count)
    except ValueError, e:
        print >> sys.stderr,  "WARNING: out of dynamic count range"
        exit(1)
    #print lower, upper
    if lower == upper:
        mean_interp = lower[1]
    elif upper[0] != lower[0]:
        slope = (upper[1] - lower[1]) / float((upper[0]-lower[0]))
        mean_interp = lower[1] + slope*(count - lower[0])
    else:
        print >> sys.stderr, "Warning: multiple empirical overdispersion values for %f" % count
        mean_interp = upper[1]
    if mean_interp < count:
        mean_interp = count
    return mean_interp
    #print lower, upper
                
# Creates a replicate, using a dispersion model
# the scale_factor scales the mean fragments requested to match
# the scale of the dispersion model (so we disperse based on a common scale of expression)
def create_replicate(original_sample_filename, rep_filename, total_frags,  disp_table, scale_factor):
    
    print >> sys.stderr, "[%s] Creating a replicate of %s in %s" % (right_now(), original_sample_filename, rep_filename)
    global overdispersion
    global poisson
    if overdispersion == 0.0 and disp_table == None:
        poisson = True
        print >> sys.stderr, "\tusing poisson dispersion model"
    else:
        print >> sys.stderr, "\tusing overdispersion model"
         
    original_sample = open(original_sample_filename, "r")
    rep = open(rep_filename, "w")
    header = original_sample.readline()
    
    print >> sys.stderr, "orig_cov\torig_frags\tnew_cov\tnew_frags"
    total_rep_frags = 0.0
    old_frags = 0.0
    total_fpkm = 0.0
    new_expr_recs = []
    for line in original_sample:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) < 7:
            print >> sys.stderr, "Warning: bad line in file, skipping"
            continue
        frag_cov = float(cols[5])
        eff_len = float(cols[6])
        fpkm = float(cols[2])
        rho = float(cols[3])
        
        num_frags = fpkm * (eff_len/1000.0) * (total_frags/1000000.0)
        old_frags += num_frags
        
        if num_frags >= 1:
            if disp_table == None:
                rep_num_frags = numpy.random.poisson(num_frags ** (overdispersion))
            else:
                # we want to maintain the index of dispersion (variance/mean) to look just like 
                # the empirically-derived model
                dispersion_frag_mean = num_frags * scale_factor
                var = get_empirical_var(disp_table, dispersion_frag_mean) / float(scale_factor)
                
                # Scale up the index of dispersion by a constant
                var *= overdispersion
                
                #var *= 1.0/(scale_factor * scale_factor) 
                #rep_num_frags = numpy.random.poisson(var)
                if var > num_frags:
                    r = (num_frags**2.0)/(var - num_frags)
                    p = num_frags / var
                    #print "num_frags = %lf, var = %lf, p = %lf, r = %lf" % (num_frags, var, p, r)
                    rep_num_frags = numpy.random.negative_binomial(r,p)
                else:
                    rep_num_frags = numpy.random.poisson(num_frags)
                #print "num frags = %lf, num frags variance = %lf, drawn = %d" % (num_frags, var, rep_num_frags)
                
        else:
            rep_num_frags = 0
        
        
        total_rep_frags += rep_num_frags
        #print >> sys.stderr, "%f\t%f" % (num_frags, rep_num_frags)
        new_expr_recs.append([rep_num_frags, cols])
    print >> sys.stderr, "Records in replicate: ", len(new_expr_recs)
    
    for rec in new_expr_recs:
        new_frags = rec[0]
        eff_len = float(rec[1][6])
        if new_frags >0 and eff_len > 0:
            new_fpkm = new_frags / (eff_len / 1000.0) / (total_rep_frags/1000000.0) 
        else:
            new_fpkm = 0
        
        total_fpkm += new_fpkm
        rec[0] = new_fpkm
        
    print >> rep, header,
    #print >> sys.stderr, new_expr_recs
    
    print >> sys.stderr, "Total original fragments: ", old_frags
    print >> sys.stderr, "Total replicate fragments: ", total_rep_frags
    for rec in new_expr_recs:
        cols = rec[1]
        eff_len = float(cols[6])
        if total_fpkm > 0:
            new_rho = rec[0] / total_fpkm
        else:
            new_rho = 0.0
            
        #print >> sys.stderr, "%f\t%f\t%f\t%f" % (old_cov, num_frags, new_frag_cov, rec[0])
        print >> rep, "%s\t%s\t%s\t%g\t%s\t%s\t%s" % (cols[0], cols[1], rec[0], new_rho, cols[4], cols[5], cols[6]) 
    rep.flush()
    rep.close()
    return rep_filename
    
    
def htseq_count(ref_gtf, sam, outfile, mode="union"):
    print >> sys.stderr, "[%s] Counting reads in %s with htseq-count" % (right_now(), sam)

    cmd = ["htseq-count"]
    cmd.append("--stranded=no")
    cmd.append("--mode=%s" % mode)
    cmd.append("-q")
    cmd.append(sam)
    cmd.append(ref_gtf)
    htseq_out = open(outfile, "w")
    
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd, stdout=htseq_out)
        if ret != 0:
            print >> sys.stderr, fail_str, "Error: could not execute htseq-count"
            exit(1)
    # tuxsim not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: htseq-count not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
    htseq_out.flush()
    htseq_out.close()

def get_fpkm_quantile(quantile_array, fpkm):
   if fpkm == 0 or fpkm < quantile_array[0]:
       return 0
   i = 1
   while i < len(quantile_array):
       if fpkm <= quantile_array[i] and fpkm > quantile_array[i-1]:
           return i
       i += 1
   return i
    
def deseq(ref_gtf, outdir, C1_replicate_sams, C2_replicate_sams, gene_fpkm_quantiles, gene_to_fpkm):
    C1_counts = {}
    C1_sams = C1_replicate_sams.split(",")
    if os.path.exists(outdir):
        pass
    else:
        os.mkdir(outdir)
    #print C1_sams
    
    global deseq_mode
    htseq_count_mode = deseq_mode
    
    for c1_sam in C1_sams:
        c1 = c1_sam.split("/")[-1]
        c1_htseqout = outdir+"/"+c1+".htseqcnt"
        htseq_count(ref_gtf, c1_sam, c1_htseqout, htseq_count_mode)
        c1_out = open(c1_htseqout)
        
        for line in c1_out:
            line = line.strip()
            cols = line.split("\t")
            gene = cols[0].strip()
            #print gene
            if gene in ["no_feature",
                        "ambiguous",
                        "too_low_aQual",
                        "not_aligned",
                        "alignment_not_unique"]:
                continue
            C1_cnts_for_gene = C1_counts.setdefault(gene, [])
            C1_cnts_for_gene.append(int(cols[1]))
            
    C2_counts = {}
    C2_sams = C2_replicate_sams.split(",")
    print C2_sams
    for c2_sam in C2_sams:
        c2 = c2_sam.split("/")[-1]
        c2_htseqout = outdir+"/"+c2+".htseqcnt"
        htseq_count(ref_gtf, c2_sam, c2_htseqout, htseq_count_mode)
        c2_out = open(c2_htseqout)
        
        for line in c2_out:
            line = line.strip()
            cols = line.split("\t")
            gene = cols[0].strip()
            #print gene
            
            if gene in ["no_feature",
                        "ambiguous",
                        "too_low_aQual",
                        "not_aligned",
                        "alignment_not_unique"]:
                continue
            C2_cnts_for_gene = C2_counts.setdefault(gene, [])
            C2_cnts_for_gene.append(int(cols[1]))
            
    header = "gene\t"
    header += "\t".join(["C1_%d" % i for i in range(0,len(C1_sams))])
    header += "\t"
    header += "\t".join(["C2_%d" % i for i in range(0,len(C2_sams))])
    
    outfile = open(outdir+"/C1_C2_htseq.counts", "w")
    
    print >> outfile, header
    for key, C1_vals in C1_counts.iteritems():
        C2_vals = C2_counts[key]
        print >> outfile, key,
        for v in C1_vals:
            print >> outfile, "\t%d" % v,
        for v in C2_vals:
            print >> outfile, "\t%d" % v,
        print >> outfile, "\n",
    outfile.flush()
    outfile.close()
    
    deseq_out = tmp_dir +"L1/perturbed_diffs/deseq/C1_C2.deseq"
    
    # Include the script inline for self-containment of the test pipeline
    R_script_str = """    
        library(DESeq)

        countsTable <- read.table("%s", header=T)
        rownames(countsTable) <- countsTable$gene
        countsTable <- countsTable[,-1]
        head(countsTable)

        conds <- c(rep("C1", ncol(countsTable)/2), rep("C2", ncol(countsTable)/2))


        cds <- newCountDataSet( countsTable, conds )
        cds <- estimateSizeFactors( cds )
        cds <- estimateDispersions( cds )
        res <- nbinomTest( cds, "C1", "C2")
        res <-res[!is.nan(res$foldChange),]
        #resSig <-res[res$padj<0.05,]
        
        write.table(res, "%s", quote=FALSE, row.names=F, sep="\\t")""" % (tmp_dir +"L1/perturbed_diffs/deseq/C1_C2_htseq.counts", deseq_out)
    
    R_script = open(tmp_dir+"deseq_script.q", "w")
    print >> R_script, R_script_str
    R_script.close()
    
    cmd = ["R", "CMD", "BATCH", "--", tmp_dir+"deseq_script.q", tmp_dir+"deseq.out"]
    try:       
        print >> run_log, " ".join(cmd)
        ret = subprocess.call(cmd)
    # R not found
    except OSError, o:
        if o.errno == errno.ENOTDIR or o.errno == errno.ENOENT:
            print >> sys.stderr, fail_str, "Error: R not found on this system.  Did you forget to include it in your PATH?"
        exit(1)
    
    significant_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    
    #print >> sys.stderr, gene_to_fpkm
    
    #sig_genes = set([])
    sig_out = open(deseq_out)
    sig_out.readline()
    for line in sig_out:
        line = line.strip()
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        gene  = cols[0].strip()
        gene_fpkm = gene_to_fpkm.get(gene)
        if gene_fpkm == None:
            #print >> sys.stderr, "DESeq: Warning, %s not found in gene fpkm table" % cols[0]
            continue
        else:
            #print >> sys.stderr, "DESeq: found %s" % cols[0]
            pass
        q_val = float(cols[7])
        if q_val <= 0.05:
            q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
            significant_genes[q].add(gene)
    print significant_genes
    return significant_genes
      
from math import modf, floor
def quantile(x, q,  qtype = 7, issorted = False):
    """
    Args:
       x - input data
       q - quantile
       qtype - algorithm
       issorted- True if x already sorted.

    Compute quantiles from input array x given q.For median,
    specify q=0.5.

    References:
       http://reference.wolfram.com/mathematica/ref/Quantile.html
       http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile

    Author:
	Ernesto P.Adorio Ph.D.
	UP Extension Program in Pampanga, Clark Field.
    """
    if not issorted:
        y = sorted(x)
    else:
        y = x
    if not (1 <= qtype <= 9): 
       return None  # error!

    # Parameters for the Hyndman and Fan algorithm
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3

            (0,   0, 0, 1), # California linear interpolation, R type 4
            (0.5, 0, 0, 1), # hydrologists method, R type 5
            (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
            (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
            (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
            (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
           ]

    a, b, c, d = abcd[qtype-1]
    n = len(x)
    g, j = modf( a + (n+b) * q -1)
    if j < 0:
        return y[0]
    elif j >= n:           
        return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!

    j = int(floor(j))
    if g ==  0:
       return y[j]
    else:
       return y[j] + (y[j+1]- y[j])* (c + d * g)
           
           
def run_cuffdiff_experiment(config, params, control, perturbed, depth_label, depth_fraction, total_frags, num_replicates, disp_table, scale_factor):

    if os.path.exists(tmp_dir + depth_label):
        pass
    else:        
        os.mkdir(tmp_dir + depth_label)
    
    ref_gtf = config.get('source_pool', 'mrna_gtf')
    
    new_config = copy.copy(config)
    num_frags = int(new_config.get('sequencing', 'num_fragments'))
    num_frags *= depth_fraction
    new_config.set('sequencing', 'num_fragments', str(int(num_frags)))
    new_config_filename = tmp_dir + depth_label + "/" + "sim.cfg" 
    new_config.write(open(new_config_filename, "w"))
    
    C1_replicate_sams = []
    C2_replicate_sams = []
    
    for i in range(1,num_replicates+1):
        # Create the initial control replicate
        C1_sam_filename = tmp_dir + depth_label + ("/C1_R%d.sam" % i)
        C1_sim_expr_filename = tmp_dir + depth_label + ("/C1_R%d_input.simexpr" % i)
        create_replicate(control, C1_sim_expr_filename, total_frags, disp_table, scale_factor)
        tuxsim(new_config_filename, params, C1_sim_expr_filename, tmp_dir + depth_label + ("/C1_R%d" % i))
        C1_replicate_sams.append(C1_sam_filename)

            
        # Create the first perturbed condition replicate (C2, L1, R1)
        C2_sam_filename = tmp_dir + depth_label + ("/C2_R%d.sam" % i)
        C2_sim_expr_filename = tmp_dir + depth_label + ("/C2_R%d_input.simexpr" % i)
        create_replicate(perturbed, C2_sim_expr_filename, total_frags, disp_table, scale_factor)
        tuxsim(new_config_filename, params, C2_sim_expr_filename, tmp_dir + depth_label + ("/C2_R%d" % i))
        C2_replicate_sams.append(C2_sam_filename)
    
    C1_replicate_sams = ",".join(C1_replicate_sams)
    C2_replicate_sams = ",".join(C2_replicate_sams)
    
    global experiment_type
    global read_length
    min_reads = min_frags_for_perturbation
    # Are we doing a read length nomogram?
    if experiment_type == "length-curve":
        # Run the differential analysis
        cuffdiff(params, 
                 tmp_dir+depth_label+"/perturbed_diffs/L1.0",
                 ref_gtf,
                 C1_replicate_sams, 
                 C2_replicate_sams,
                 ["C1", "C2"],
                 #0,
                 min_reads,
                 read_skip_fraction=0.5,
                 no_read_pairing=True,
                 read_trim_length=floor(1.0 * read_length),
                 library_type="transfrags" 
                 )
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/L0.75",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.5,
                no_read_pairing=True,
                read_trim_length=floor(0.75 * read_length),
                 library_type="transfrags"
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/L0.50",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.5,
                no_read_pairing=True,
                read_trim_length=floor(0.5 * read_length),
                 library_type="transfrags"
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/L0.25",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.5,
                no_read_pairing=True,
                read_trim_length=floor(0.25 * read_length),
                library_type="transfrags"
                )
            
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/L0.15",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.5,
                no_read_pairing=True,
                read_trim_length=floor(0.15 * read_length),
                library_type="transfrags"
                )
            
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/L0.05",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.5,
                no_read_pairing=True,
                read_trim_length=floor(0.05 * read_length),
                library_type="transfrags"
                )
    elif experiment_type == "depth-curve":   
        # Run the differential analysis for the depth of sequencing nomogram
        # Do the paired samples
        cuffdiff(params, 
                 tmp_dir+depth_label+"/perturbed_diffs/r1.0",
                 ref_gtf,
                 C1_replicate_sams, 
                 C2_replicate_sams,
                 ["C1", "C2"],
                 #0,
                 min_reads,
                 read_skip_fraction=0.0,
                 no_read_pairing=False
                 )
              
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.75",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.25,
                no_read_pairing=False
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.50",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.50,
                no_read_pairing=False
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.25",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.75,
                no_read_pairing=False
                )
            
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.15",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.85,
                no_read_pairing=False
                )
            
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.05",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.95,
                no_read_pairing=False
                )
            
        # Do the unpaired samples
        cuffdiff(params, 
                  tmp_dir+depth_label+"/perturbed_diffs/r1.0_no_pair",
                  ref_gtf,
                  C1_replicate_sams, 
                  C2_replicate_sams,
                  ["C1", "C2"],
                  #0,
                  min_reads,
                  read_skip_fraction=0.0,
                  no_read_pairing=True
                  )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.75_no_pair",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.50,
                no_read_pairing=True
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.50_no_pair",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.75,
                no_read_pairing=True
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.25_no_pair",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.875,
                no_read_pairing=True
                )

        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/r0.15_no_pair",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.925,
                no_read_pairing=True
                )

        cuffdiff(params, 
                 tmp_dir+depth_label+"/perturbed_diffs/r0.05_no_pair",
                 ref_gtf,
                 C1_replicate_sams, 
                 C2_replicate_sams,
                 ["C1", "C2"],
                 #0,
                 min_reads,
                 read_skip_fraction=0.975,
                 no_read_pairing=True
                 )
    elif experiment_type == "deseq":
        cuffdiff(params, 
                 tmp_dir+depth_label+"/perturbed_diffs/r1.0",
                 ref_gtf,
                 C1_replicate_sams, 
                 C2_replicate_sams,
                 ["C1", "C2"],
                 #0,
                 min_reads,
                 read_skip_fraction=0.0,
                 no_read_pairing=False
                 )
    elif experiment_type == "compress-gtf":
        
        # Generate the compressed representative GTF files:
        union_gtf = output_dir+"union.gtf"
        intersection_gtf = output_dir+"intersection.gtf"
        
        compress_gtf(params, ref_gtf, union_gtf, "--union")
        compress_gtf(params, ref_gtf, intersection_gtf, "--intersection")
        
        # Run Cuffdiff in all three modes for downstream comparison
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/normal",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.0,
                no_read_pairing=False
                )
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/union",
                union_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.0,
                no_read_pairing=False
                )
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/intersection",
                intersection_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.0,
                no_read_pairing=False
                )
    
    elif experiment_type == "basic":        
        # Run Cuffdiff in all normal modes for downstream comparison
        cuffdiff(params, 
                tmp_dir+depth_label+"/perturbed_diffs/normal",
                ref_gtf,
                C1_replicate_sams, 
                C2_replicate_sams,
                ["C1", "C2"],
                #0,
                min_reads,
                read_skip_fraction=0.0,
                no_read_pairing=False
                )
                
    elif experiment_type == "replicate-series":        
        C1_reps = C1_replicate_sams.split(",")
        C2_reps = C2_replicate_sams.split(",")
        if len(C1_reps) != len(C2_reps):
            print >> sys.stderr, "Error: must have equal # of replicates for both conditions"
            exit(1)
        for i in range(0, len(C1_reps)):
            C1_rep_list = C1_reps[0:i + 1]
            C2_rep_list = C2_reps[0:i + 1]
            
            #skip_fraction = float(i)/float(len(C1_reps))
            skip_fraction = 0.0
            cuffdiff(params, 
                    tmp_dir+depth_label+"/perturbed_diffs/with_%d_reps" % (i+1),
                    ref_gtf,
                    ",".join(C1_rep_list), 
                    ",".join(C2_rep_list),
                    ["C1", "C2"],
                    #0,
                    min_reads,
                    read_skip_fraction=skip_fraction,
                    no_read_pairing=False
                    )
    else:
        print >> sys.stderr, "Error: unrecognized experiment type"

    all_expressed_file = open(control)
    all_expressed_file.readline()

    iso_fpkms = []
    gene_to_fpkm = {}
    
    
    for line in all_expressed_file:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) > 3:
            fpkm = float(cols[2])
            if fpkm >= 0:
                iso_fpkms.append(fpkm)
            fpkms = gene_to_fpkm.setdefault(cols[0], [])
            fpkms.append(fpkm)
            gene_to_fpkm[cols[0]] = fpkms
    print >> sys.stderr, "HAVE %d expressed genes" % len(gene_to_fpkm)
    #print >> sys.stderr, gene_to_fpkm            
    gene_fpkms = []
    for gene, fpkms in gene_to_fpkm.iteritems():
        gene_fpkm = sum(fpkms)
        gene_to_fpkm[gene] = gene_fpkm
        if gene_fpkm > 0:
            gene_fpkms.append(gene_fpkm)
    
        
    gene_fpkm_quantiles = [quantile(gene_fpkms, x) for x in [0.25, 0.5, 0.75, 1.0]]
    

    def get_sig(cuff_out_dir, gene_fpkm_quantiles, gene_to_fpkm):
                
        significant_transcripts = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
        significant_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
        significant_spliced_tss_groups = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
        significant_promoter_switch_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
        
        isoform_exp = open(cuff_out_dir + "/isoform_exp.diff")
        isoform_exp.readline()
        for line in isoform_exp:
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 13:
                continue
            gene_fpkm = gene_to_fpkm.get(cols[1])
            if gene_fpkm == None:
                #print >> sys.stderr, "Warning %s not found in gene fpkm table" % cols[1]
                continue
            if cols[13] == "yes":
                q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
                significant_transcripts[q].add(cols[0])
                
        gene_exp = open(cuff_out_dir + "/gene_exp.diff")
        gene_exp.readline()
        for line in gene_exp:
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 13:
                continue
            gene_fpkm = gene_to_fpkm.get(cols[1])
            if gene_fpkm == None:
                #print >> sys.stderr, "Warning %s not found in gene fpkm table" % cols[1]
                continue
            if cols[13] == "yes":
                q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
                significant_genes[q].add(cols[0])
                
        splicing = open(cuff_out_dir + "/splicing.diff")
        splicing.readline()
        for line in splicing:
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 13:
                continue
            gene_fpkm = gene_to_fpkm.get(cols[1])
            if gene_fpkm == None:
                #print >> sys.stderr, "Warning %s not found in gene fpkm table" % cols[1]
                continue
            if cols[13] == "yes":
                q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
                significant_spliced_tss_groups[q].add(cols[0])
                
        promoters = open(cuff_out_dir + "/promoters.diff")
        promoters.readline()
        for line in promoters:
            line = line.strip()
            cols = line.split('\t')
            if len(cols) < 13:
                continue
            gene_fpkm = gene_to_fpkm.get(cols[1])
            if gene_fpkm == None:
                #print >> sys.stderr, "Warning %s not found in gene fpkm table" % cols[1]
                continue
            if cols[13] == "yes":
                q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
                significant_promoter_switch_genes[q].add(cols[0])
        
        return (significant_genes, significant_transcripts, significant_spliced_tss_groups, significant_promoter_switch_genes)
    
    union_sig_trans = set([])
    intersection_sig_trans = set([])
    
    global diff_mode
    perturbed_lines = open(output_dir+"/perturbed.lines")
    perturbed_lines.readline()
    true_perturbed_trans = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    true_perturbed_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    true_diff_spliced_tss_groups = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    true_promoter_switch_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    
    tss_ids_for_perturbed_genes = {}
    isoforms_for_perturbed_tss_groups = {}
    
    for line in perturbed_lines:
        line = line.strip()
        cols = line.split('\t')
        #print >> sys.stderr, len(cols), cols
        if len(cols) >0:
            gene_fpkm = gene_to_fpkm.get(cols[0])
            if gene_fpkm == None:
                print >> sys.stderr, "Warning: no gene fpkm entry for %s" % cols[0]
                continue
            q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
            if diff_mode not in  ["rotate-rhos", "reflect-rhos"]:
                true_perturbed_genes[q].add(cols[0])
            true_perturbed_trans[q].add(cols[1])
            #print len(cols)
            #print cols
            if len(cols) >= 8:
                tss_id = cols[7]
                #print tss_id
                if tss_id != "-":
                    #true_diff_spliced_tss_groups.add(tss_id)
                    isos = isoforms_for_perturbed_tss_groups.setdefault(tss_id, set([]))
                    isos.add(cols[1])
                    isoforms_for_perturbed_tss_groups[tss_id] = isos
                    
                    tss_ids = tss_ids_for_perturbed_genes.setdefault(cols[0], set([]))
                    tss_ids.add(tss_id)
                    tss_ids_for_perturbed_genes[cols[0]] = tss_ids
    
    all_expressed_file = open(control)
    all_expressed_file.readline()
    all_expressed_trans = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    all_expressed_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    all_expressed_tss_ids = [set([]) for x in range(0, len(gene_fpkm_quantiles))]

    tssid_to_gene = {}

    for line in all_expressed_file:
        line = line.strip()
        cols = line.split('\t')
        if len(cols) > 3 and float(cols[2]) > 0:
            gene_fpkm = gene_to_fpkm.get(cols[0])
            if gene_fpkm == None:
                continue
            q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
            if q >= len(gene_fpkm_quantiles):
                print >> sys.stderr, "Error: %d is out of range %s, %f"  % (q, str(gene_fpkm_quantiles), gene_fpkm)
            all_expressed_genes[q].add(cols[0])
            all_expressed_trans[q].add(cols[1])

            if len(cols) >= 8:
                tss_id = cols[7]
                if tss_id not in "-":
                    tssid_to_gene[tss_id] = cols[0]
                    all_expressed_tss_ids[q].add(tss_id)
                    isos = isoforms_for_perturbed_tss_groups.get(tss_id)
                    if isos != None:
                        isos.add(cols[1])
                        isoforms_for_perturbed_tss_groups[tss_id] = isos

                    tss_ids = tss_ids_for_perturbed_genes.get(cols[0])
                    if tss_ids != None:
                        tss_ids.add(tss_id)
                        tss_ids_for_perturbed_genes[cols[0]] = tss_ids

    for g_id, tss_ids in tss_ids_for_perturbed_genes.iteritems():
        if len(tss_ids) > 1:
            gene_fpkm = gene_to_fpkm.get(g_id)
            if gene_fpkm == None:
                continue
            q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
            true_promoter_switch_genes[q].add(g_id) 

    for tss_id, isos in isoforms_for_perturbed_tss_groups.iteritems():
        if len(isos) > 1:
            g_id = tssid_to_gene.get(tss_id)
            if g_id == None:
                print >> sys.stderr, "Warning: no gene_id for tss_id %s" % tss_id
                continue
            gene_fpkm = gene_to_fpkm.get(g_id)
            if gene_fpkm == None:
                continue
            q = get_fpkm_quantile(gene_fpkm_quantiles, gene_fpkm)
            true_diff_spliced_tss_groups[q].add(tss_id)
    
    # Initialize the sets of genes we didn't perturb
    not_perturbed_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    not_perturbed_trans = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    not_spliced_tss_groups = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    not_promoter_switch_genes = [set([]) for x in range(0, len(gene_fpkm_quantiles))]
    for i in range(0,len(gene_fpkm_quantiles)):
        not_perturbed_genes[i] = all_expressed_genes[i] - true_perturbed_genes[i]
        not_perturbed_trans[i] = all_expressed_trans[i] - true_perturbed_trans[i]
        not_spliced_tss_groups[i] = all_expressed_tss_ids[i] - true_diff_spliced_tss_groups[i]
        not_promoter_switch_genes[i] = all_expressed_genes[i] - true_promoter_switch_genes[i]
     
    #print >> sys.stderr,  true_perturbed_trans
    # Are we doing a read length nomogram?
    if experiment_type == "length-curve":
        print >> sys.stderr, "Processing results from Read Length nomogram"
        # Paired tests:
        L1_0_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L1.0/", gene_fpkm_quantiles, gene_to_fpkm)
        L0_75_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.75/", gene_fpkm_quantiles, gene_to_fpkm)
        L0_50_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.50/", gene_fpkm_quantiles, gene_to_fpkm)
        L0_25_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.25/", gene_fpkm_quantiles, gene_to_fpkm)
        L0_15_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.15/", gene_fpkm_quantiles, gene_to_fpkm)
        L0_05_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.05/", gene_fpkm_quantiles, gene_to_fpkm)

        # # Unpaired tests:
        # L1_0_no_pair_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L1.0_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        # L0_75_no_pair_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.75_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        # L0_50_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.50_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        # L0_25_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.25_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        # L0_15_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.15_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        # L0_05_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/L0.05_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
    
        exps =  [("normal", False, 1.0, 1.00, num_replicates, L1_0_acc[0], L1_0_acc[1], L1_0_acc[2], L1_0_acc[3]),
                 ("normal", False, 1.0, 0.75, num_replicates, L0_75_acc[0], L0_75_acc[1], L0_75_acc[2], L0_75_acc[3]),
                 ("normal", False, 1.0, 0.50, num_replicates, L0_50_acc[0], L0_50_acc[1], L0_50_acc[2], L0_50_acc[3]),
                 ("normal", False, 1.0, 0.25, num_replicates, L0_25_acc[0], L0_25_acc[1], L0_25_acc[2], L0_25_acc[3]),
                 ("normal", False, 1.0, 0.15, num_replicates, L0_15_acc[0], L0_15_acc[1], L0_15_acc[2], L0_15_acc[3]),
                 ("normal", False, 1.0, 0.05, num_replicates, L0_05_acc[0], L0_05_acc[1], L0_05_acc[2], L0_05_acc[3])]
                 
    elif experiment_type == "depth-curve": # Then it's a depth of sequencing curve 
        print >> sys.stderr, "Processing results from Sequencing depth nomogram"   
        # Paired tests:
        r1_0_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r1.0/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_75_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.75/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_50_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.50/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_25_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.25/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_15_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.15/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_05_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.05/", gene_fpkm_quantiles, gene_to_fpkm)

        # Unpaired tests:
        r1_0_no_pair_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r1.0_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_75_no_pair_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.75_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_50_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.50_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_25_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.25_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_15_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.15_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
        r0_05_no_pair_acc  = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r0.05_no_pair/", gene_fpkm_quantiles, gene_to_fpkm)
    
        exps =  [("normal", True, 1.0, 1.0, num_replicates, r1_0_acc[0], r1_0_acc[1], r1_0_acc[2], r1_0_acc[3]),
                 ("normal", True, 0.75,1.0,  num_replicates, r0_75_acc[0], r0_75_acc[1], r0_75_acc[2], r0_75_acc[3]),
                 ("normal", True, 0.50,1.0,  num_replicates, r0_50_acc[0], r0_50_acc[1], r0_50_acc[2], r0_50_acc[3]),
                 ("normal", True, 0.25,1.0,  num_replicates, r0_25_acc[0], r0_25_acc[1], r0_25_acc[2], r0_25_acc[3]),
                 ("normal", True, 0.15,1.0,  num_replicates, r0_15_acc[0], r0_15_acc[1], r0_15_acc[2], r0_15_acc[3]),
                 ("normal", True, 0.05,1.0,  num_replicates, r0_05_acc[0], r0_05_acc[1], r0_05_acc[2], r0_05_acc[3]),
                 ("normal", False, 1.0, 1.0,  num_replicates, r1_0_no_pair_acc[0], r1_0_no_pair_acc[1], r1_0_no_pair_acc[2], r1_0_no_pair_acc[3]),
                 ("normal", False, 0.75, 1.0,  num_replicates, r0_75_no_pair_acc[0], r0_75_no_pair_acc[1], r0_75_no_pair_acc[2], r0_75_no_pair_acc[3]),
                 ("normal", False, 0.50, 1.0,  num_replicates, r0_50_no_pair_acc[0], r0_50_no_pair_acc[1], r0_50_no_pair_acc[2], r0_50_no_pair_acc[3]),
                 ("normal", False, 0.25, 1.0,  num_replicates, r0_25_no_pair_acc[0], r0_25_no_pair_acc[1], r0_25_no_pair_acc[2], r0_25_no_pair_acc[3]),
                 ("normal", False, 0.15, 1.0,  num_replicates, r0_15_no_pair_acc[0], r0_15_no_pair_acc[1], r0_15_no_pair_acc[2], r0_15_no_pair_acc[3]),
                 ("normal", False, 0.05, 1.0,  num_replicates, r0_05_no_pair_acc[0], r0_05_no_pair_acc[1], r0_05_no_pair_acc[2], r0_05_no_pair_acc[3])]
    elif experiment_type == "deseq":
        print >> sys.stderr, "Running analysis with DESeq"   
        deseq_sig_genes = deseq(ref_gtf, 
                                tmp_dir+depth_label+"/perturbed_diffs/deseq", 
                                C1_replicate_sams, 
                                C2_replicate_sams, 
                                gene_fpkm_quantiles, 
                                gene_to_fpkm)
        #print deseq_sig_genes
        deseq_acc = [deseq_sig_genes, 
                     [set([]) for x in range(0, len(gene_fpkm_quantiles))],
                     [set([]) for x in range(0, len(gene_fpkm_quantiles))],
                     [set([]) for x in range(0, len(gene_fpkm_quantiles))]]

        r1_0_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/r1.0/", gene_fpkm_quantiles, gene_to_fpkm)
    
        exps = [("normal", True, 1.0, 1.0, num_replicates, r1_0_acc[0], r1_0_acc[1], r1_0_acc[2], r1_0_acc[3]),
                ("deseq", True, 1.0, 1.0, num_replicates, deseq_acc[0], deseq_acc[1], deseq_acc[2], deseq_acc[3])]
                
    elif experiment_type == "compress-gtf":
         #print >> sys.stderr, "Running analysis with Compressed GTFs"   

         normal_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/normal/", gene_fpkm_quantiles, gene_to_fpkm)
         union_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/union/", gene_fpkm_quantiles, gene_to_fpkm)
         intersection_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/intersection/", gene_fpkm_quantiles, gene_to_fpkm)

         exps =  [("normal", True, 1.0, 1.0, num_replicates, r1_0_acc[0], r1_0_acc[1], r1_0_acc[2], r1_0_acc[3]),
                  ("union", True, 1.0, 1.0, num_replicates, union_acc[0], union_acc[1], union_acc[2], union_acc[3]),
                  ("intersection", True, 1.0, 1.0, num_replicates, intersection_acc[0], intersection_acc[1], intersection_acc[2], intersection_acc[3])]
    elif experiment_type == "basic":
        print >> sys.stderr, "Running basic analysis"   
        normal_acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/normal/", gene_fpkm_quantiles, gene_to_fpkm)
        exps =  [("normal", True, 1.0, 1.0, num_replicates, normal_acc[0], normal_acc[1], normal_acc[2], normal_acc[3])]
    elif experiment_type == "replicate-series":        
        C1_reps = C1_replicate_sams.split(",")
        C2_reps = C2_replicate_sams.split(",")
        if len(C1_reps) != len(C2_reps):
            print >> sys.stderr, "Error: must have equal # of replicates for both conditions"
            exit(1)
        exps =  []
        for i in range(0, len(C1_reps)):
            acc = get_sig(tmp_dir+depth_label+"/perturbed_diffs/with_%d_reps/" % (i+1), gene_fpkm_quantiles, gene_to_fpkm)
            exps.append(("normal", True, 1.0, 1.0, i + 1, acc[0], acc[1], acc[2], acc[3]))

    
    def get_accuracy(sig_genes, sig_trans, sig_splice, sig_prom,
                     true_perturbed_genes, true_perturbed_trans, true_diff_spliced_tss_groups, true_promoter_switch_genes,
                     not_perturbed_genes, not_perturbed_trans, not_spliced_tss_groups, not_promoter_switch_genes):
        g_tp = sig_genes & true_perturbed_genes
        g_tn = not_perturbed_genes - sig_genes
        g_fp = sig_genes & not_perturbed_genes
        g_fn = true_perturbed_genes - sig_genes
    
        t_tp = sig_trans & true_perturbed_trans
        t_tn = not_perturbed_trans - sig_trans
        t_fp = sig_trans & not_perturbed_trans
        t_fn = true_perturbed_trans - sig_trans
    
        s_tp = sig_splice & true_diff_spliced_tss_groups
        s_tn = not_spliced_tss_groups - sig_splice
        s_fp = sig_splice & not_spliced_tss_groups
        s_fn = true_diff_spliced_tss_groups - sig_splice
    
        p_tp = sig_prom & true_promoter_switch_genes
        p_tn = not_promoter_switch_genes - sig_prom
        p_fp = sig_prom & not_promoter_switch_genes
        p_fn = true_promoter_switch_genes - sig_prom

        if (len(g_tp) + len(g_fp)) > 0:
            g_prec = len(g_tp) / float(len(g_tp) + len(g_fp))
        else:
            g_prec = 0.0
        if (len(g_tp) + len(g_fn) > 0):
            g_recall = len(g_tp) / float(len(g_tp) + len(g_fn))
        else:
            g_recall = 0.0
    
        if (len(t_tp) + len(t_fp)) > 0:
            t_prec = len(t_tp) / float(len(t_tp) + len(t_fp))
        else:
            t_prec = 0.0
        if (len(t_tp) + len(t_fn) > 0):
            t_recall = len(t_tp) / float(len(t_tp) + len(t_fn))
        else:
            t_recall = 0.0
        
        if (len(s_tp) + len(s_fp)) > 0:
            s_prec = len(s_tp) / float(len(s_tp) + len(s_fp))
        else:
            s_prec = 0.0
        if (len(s_tp) + len(s_fn) > 0):
            s_recall = len(s_tp) / float(len(s_tp) + len(s_fn))
        else:
            s_recall = 0.0
        
        if (len(p_tp) + len(p_fp)) > 0:
            p_prec = len(p_tp) / float(len(p_tp) + len(p_fp))
        else:
            p_prec = 0.0
        if (len(p_tp) + len(p_fn) > 0):
            p_recall = len(p_tp) / float(len(p_tp) + len(p_fn))
        else:
            p_recall = 0.0
        return (g_prec, g_recall, t_prec, t_recall, s_prec, s_recall, p_prec, p_recall)
    
    diff_acc = open(output_dir+"diff_accuracy.out", "w")
    print >> diff_acc, "Mode\tFPKM_quantile\tPaired\tExperiment\tDepth\tReadLength\tNumReplicates\tGenePrec\tGeneRecall\tTransPrec\tTransRecall\tSplicingPrec\tSplicingRecall\tPromotersPrec\tPromotersRecall"
    for (mode, paired, depth, read_length, n_reps, sig_genes, sig_trans, sig_splice, sig_prom) in exps:
        for i in range(0, len(sig_genes)):
            accuracy = get_accuracy(sig_genes[i], sig_trans[i], sig_splice[i], sig_prom[i],
                                    true_perturbed_genes[i], true_perturbed_trans[i], true_diff_spliced_tss_groups[i], true_promoter_switch_genes[i],
                                    not_perturbed_genes[i], not_perturbed_trans[i], not_spliced_tss_groups[i], not_promoter_switch_genes[i])
            (g_prec, g_recall, t_prec, t_recall, s_prec, s_recall, p_prec, p_recall) = accuracy
            print >> diff_acc, "%s\t%d\t%s\t%s\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (mode, i, str(paired), diff_mode, depth, read_length, n_reps, g_prec, g_recall, t_prec, t_recall, s_prec, s_recall, p_prec, p_recall)   
        all_sig_genes  = set([item for sublist in sig_genes for item in sublist])
        all_sig_trans  = set([item for sublist in sig_trans for item in sublist])
        all_sig_splice = set([item for sublist in sig_splice for item in sublist])
        all_sig_prom   = set([item for sublist in sig_prom for item in sublist])
        all_true_sig_genes  = set([item for sublist in true_perturbed_genes for item in sublist])
        all_true_sig_trans  = set([item for sublist in true_perturbed_trans for item in sublist])
        all_true_sig_splice = set([item for sublist in true_diff_spliced_tss_groups for item in sublist])
        all_true_sig_prom   = set([item for sublist in true_promoter_switch_genes for item in sublist])
        all_not_sig_genes  = set([item for sublist in not_perturbed_genes for item in sublist])
        all_not_sig_trans  = set([item for sublist in not_perturbed_trans for item in sublist])
        all_not_sig_splice = set([item for sublist in not_spliced_tss_groups for item in sublist])
        all_not_sig_prom   = set([item for sublist in not_promoter_switch_genes for item in sublist])
        accuracy = get_accuracy(all_sig_genes, all_sig_trans, all_sig_splice, all_sig_prom,
                                all_true_sig_genes, all_true_sig_trans, all_true_sig_splice, all_true_sig_prom,
                                all_not_sig_genes, all_not_sig_trans, all_not_sig_splice, all_not_sig_prom)
        (g_prec, g_recall, t_prec, t_recall, s_prec, s_recall, p_prec, p_recall) = accuracy
        print >> diff_acc, "%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" % (mode, "all", str(paired), diff_mode, depth, read_length, n_reps, g_prec, g_recall, t_prec, t_recall, s_prec, s_recall, p_prec, p_recall)
    diff_acc.flush()
    diff_acc.close()
    
def main(argv=None):
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
    
    # Initialize default parameter values
    params = TestParams()
    
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()
            
        print >> sys.stderr
        print >> sys.stderr, "[%s] Beginning Multiplexed RNA-Seq test (suite v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------" 
        
        global iso_switch
        
        start_time = datetime.now()
        prepare_output_dir()
        
        global run_log
        run_log = open(logging_dir + "run.log", "w", 0)
        global run_cmd
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd
        tuxsim_cfg_params = {}

        config = ConfigParser.SafeConfigParser()
        print >> sys.stderr, args[0]
        config.read(args[0])
        
        ref_gtf = config.get('source_pool', 'mrna_gtf')
        total_frags = float(config.get('sequencing', 'num_fragments'))
        
        global read_length
        read_length = float(config.get('sequencing', 'read_length'))
        global frag_length_mean
        frag_length_mean = float(config.get('fragment', 'length.mean'))
        global frag_length_std_dev
        frag_length_std_dev = float(config.get('fragment', 'length.std_dev'))
        
        disp_table = None
        scale_factor = 1.0
        if dispersion_filename != None:
            (disp_table, frags_in_table) = load_empirical_dispersion(dispersion_filename)
            scale_factor = frags_in_table / total_frags
        
        # config.set('sequencing', 'num_fragments', "1")
        # config.write(open('sim2.cfg', "w"))
        # return 
        
        config_file = args[0]
        # Set up the initial control condition at full depth (C1, L1, R1)
        #base_sam_filename = "sim.sam"
        if extern_simexpr_filename != None:
            base_sim_expr_filename = extern_simexpr_filename    
            #tuxsim("sim.cfg", params, extern_simexpr_filename)  
        else:  
            base_sim_expr_filename = output_dir + "sim.simexpr"    
            tuxsim(args[0], params, out_prefix=output_dir+"sim")
        
        global min_fold_change
        global max_fold_change
        
        # Set up the perturbed condition at full depth 
        #perturbed_sam_filename = "perturbed.sam"
        
        global diff_mode
        #print diff_mode
        
        global num_perturbed_genes
        
        perturbed_sim_expr_filename = output_dir+"sim_perturbed.simexpr"      
        if diff_mode == "iso-switch":
            create_isoform_switch_condition(base_sim_expr_filename, 
                                            perturbed_sim_expr_filename, 
                                            output_dir+"/perturbed.lines", 
                                            max_fold_change, 
                                            num_perturbed_genes, 
                                            total_frags)
        elif diff_mode == "single-iso-up":
            create_single_iso_up_condition(base_sim_expr_filename, 
                                           perturbed_sim_expr_filename, 
                                           output_dir+"/perturbed.lines", 
                                           max_fold_change, 
                                           num_perturbed_genes, 
                                           total_frags)
        elif diff_mode == "interval":
            create_interval_fold_change_condition(base_sim_expr_filename, 
                                          perturbed_sim_expr_filename, 
                                          output_dir+"/perturbed.lines", 
                                          1,
                                          max_fold_change, 
                                          num_perturbed_genes, 
                                          total_frags,
                                          "single")
        elif diff_mode == "interval-multi":
            create_interval_fold_change_condition(base_sim_expr_filename, 
                                        perturbed_sim_expr_filename, 
                                        output_dir+"/perturbed.lines", 
                                        min_fold_change,
                                        max_fold_change, 
                                        num_perturbed_genes, 
                                        total_frags,
                                        "multi")
        elif diff_mode == "uniform":
            create_perturbed_condition(base_sim_expr_filename, 
                                       perturbed_sim_expr_filename, 
                                       output_dir+"/perturbed.lines", 
                                       max_fold_change, 
                                       num_perturbed_genes, 
                                       total_frags)
        elif diff_mode == "shift-to-single":
            create_shift_all_to_one_condition(base_sim_expr_filename, 
                                              perturbed_sim_expr_filename, 
                                              output_dir+"/perturbed.lines", 
                                              max_fold_change, 
                                              num_perturbed_genes, 
                                              total_frags)
        elif diff_mode == "rotate-rhos":
            create_rotate_rhos_condition(base_sim_expr_filename, 
                                         perturbed_sim_expr_filename, 
                                         output_dir+"/perturbed.lines",  
                                         num_perturbed_genes, 
                                         total_frags)
        elif diff_mode == "reflect-rhos":
            create_reflect_isoforms_condition(base_sim_expr_filename, 
                                              perturbed_sim_expr_filename, 
                                              output_dir+"/perturbed.lines",  
                                              num_perturbed_genes, 
                                              total_frags)
        else:
            print >> sys.stderr, "Error: unrecognized perturbation mode %s" % (diff_mode)
            exit(1)
            
        #tuxsim(args[0], params, perturbed_sim_expr_filename, output_dir+"perturbed")
        
        L1_sig = run_cuffdiff_experiment(config, params, base_sim_expr_filename, perturbed_sim_expr_filename, "L1", 1.0, total_frags, num_replicates, disp_table, scale_factor)
        
        finish_time = datetime.now()
        duration = finish_time - start_time
        run_log.flush()
        run_log.close()
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Test complete [%s elapsed]" %  formatTD(duration)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        #print >> sys.stderr, "    for detailed help see http://spats.cbcb.umd.edu/manual.html"
        return 2


if __name__ == "__main__":
    sys.exit(main())
