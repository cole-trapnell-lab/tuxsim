include $(top_srcdir)/config/tuxsim.mk

AM_CPPFLAGS = -I$(top_srcdir)/src 

EXTRA_DIST = $(top_srcdir)/AUTHORS $(top_srcdir)/make_bin.sh

AM_CXXFLAGS =

bin_PROGRAMS = \
	tuxsim

noinst_HEADERS = \
	abundances.h \
	common.h \
    codons.h \
    fragments.h \
    GArgs.h \
    GBase.h \
    gdna.h \
    GFaSeqGet.h \
    GFastaFile.h \
    gff.h \
    GHash.hh \
    GStr.h \
    GList.hh \
    gtf_tracking.h \
    hits.h \
    options.h \
    scaffolds.h \
    sequencing.h 

noinst_LIBRARIES = libtuxsim.a

libtuxsim_a_SOURCES = \
	abundances.cpp \
	common.cpp \
    codons.cpp \
    fragments.cpp \
    GArgs.cpp \
    GBase.cpp \
    gdna.cpp \
    GFaSeqGet.cpp \
    gff.cpp \
    GStr.cpp \
    gtf_tracking.cpp \
    hits.cpp \
    options.cpp \
    scaffolds.cpp \
    sequencing.cpp
    
    
tuxsim_SOURCES = tuxsim.cpp  
tuxsim_LDADD = libtuxsim.a $(BOOST_THREAD_LIB) $(BOOST_PROGRAM_OPTIONS_LIB)
tuxsim_LDFLAGS =  $(BOOST_LDFLAGS)