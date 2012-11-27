/*
 *  TestHits.cpp
 *  tuxsim
 *
 *  Created by Cole Trapnell on 5/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include "TestHits.h"

#include "hits.h"

#include <cppunit/extensions/HelperMacros.h>

CPPUNIT_TEST_SUITE_REGISTRATION(TestHits);

//CppUnit::Test *TestHits::suite()
//{
//  CppUnit::TestSuite *suiteOfTests = new CppUnit::TestSuite("TestHits");
//  suiteOfTests->addTest(new CppUnit::TestCaller<TestWchar>("testMateHitLessThan",&TestWchar::testMateHitLessThan));
//  return suiteOfTests;
//}

void TestHits::testMateHitLessThan()
{

}

void TestHits::testMateHitEquivalent()
{
    
}

void TestHits::testSingletonLessThanCigarOpLen()
{
    vector<CigarOp> read_1_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_1_1(new ReadHit(1, 1, 1, read_1_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m1(1, read_1_1,  shared_ptr<ReadHit>(), 0, 0);
    
    vector<CigarOp> read_2_1_cig(1, CigarOp(MATCH, 75));
    shared_ptr<ReadHit> read_2_1(new ReadHit(1, 1, 1, read_2_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m2(1, read_2_1,  shared_ptr<ReadHit>(), 0, 0);
    
    CPPUNIT_ASSERT(m1 < m2);
    CPPUNIT_ASSERT(!(m2 < m1));
}

void TestHits::testSingletonLessThanCigarOpCode()
{
    vector<CigarOp> read_1_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_1_1(new ReadHit(1, 1, 1, read_1_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m1(1, read_1_1,  shared_ptr<ReadHit>(), 0, 0);
    
    vector<CigarOp> read_2_1_cig(1, CigarOp(REF_SKIP, 76));
    shared_ptr<ReadHit> read_2_1(new ReadHit(1, 1, 1, read_2_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m2(1, read_2_1,  shared_ptr<ReadHit>(), 0, 0);
    
    CPPUNIT_ASSERT(m1 < m2);
    CPPUNIT_ASSERT(!(m2 < m1));
}

void TestHits::testSingletonLessThanCigarLength()
{
    vector<CigarOp> read_1_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_1_1(new ReadHit(1, 1, 1, read_1_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m1(1, read_1_1,  shared_ptr<ReadHit>(), 0, 0);
    
    vector<CigarOp> read_2_1_cig;
    read_2_1_cig.push_back(CigarOp(MATCH,36));
    read_2_1_cig.push_back(CigarOp(REF_SKIP,10));
    read_2_1_cig.push_back(CigarOp(MATCH,30));
    
    shared_ptr<ReadHit> read_2_1(new ReadHit(1, 1, 1, read_2_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m2(1, read_2_1,  shared_ptr<ReadHit>(), 0, 0);
    
    CPPUNIT_ASSERT(m1 < m2);
    CPPUNIT_ASSERT(!(m2 < m1));
}

void TestHits::testSingletonLessThanInsertID()
{
    vector<CigarOp> read_1_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_1_1(new ReadHit(1, 1, 1, read_1_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m1(1, read_1_1,  shared_ptr<ReadHit>(), 0, 0);
    
    vector<CigarOp> read_2_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_2_1(new ReadHit(1, 2, 1, read_2_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m2(1, read_2_1,  shared_ptr<ReadHit>(), 0, 0);
    
    CPPUNIT_ASSERT(m1 < m2);
    CPPUNIT_ASSERT(!(m2 < m1));
}


void TestHits::testSingletonEquivalant()
{
    vector<CigarOp> read_1_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_1_1(new ReadHit(1, 1, 1, read_1_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m1(1, read_1_1,  shared_ptr<ReadHit>(), 0, 0);
    
    vector<CigarOp> read_2_1_cig(1, CigarOp(MATCH, 76));
    shared_ptr<ReadHit> read_2_1(new ReadHit(1, 1, 1, read_2_1_cig, false, CUFF_STRAND_UNKNOWN, 1, 100, 0, 0));
    
    MateHit m2(1, read_2_1,  shared_ptr<ReadHit>(), 0, 0);
    
    CPPUNIT_ASSERT(!(m1 < m2));
    CPPUNIT_ASSERT(!(m2 < m1));
}