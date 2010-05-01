/*
 *  TestHits.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 5/1/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#ifndef TESTMATEHITS_H
#define TESTMATEHITS_H

#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestSuite.h>
#include <cppunit/extensions/HelperMacros.h>

class TestHits : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(TestHits);
    //CPPUNIT_TEST(testMateHitLessThan);
    //CPPUNIT_TEST(testMateHitLessThan);
    //CPPUNIT_TEST(testMateHitLessThan);
    CPPUNIT_TEST(testSingletonLessThanCigarOpLen);
    CPPUNIT_TEST(testSingletonLessThanCigarOpCode);
    CPPUNIT_TEST(testSingletonLessThanCigarLength);
    CPPUNIT_TEST(testSingletonLessThanInsertID);
    CPPUNIT_TEST(testSingletonEquivalant);
    CPPUNIT_TEST_SUITE_END();
    
public:
    void testMateHitLessThan();
    void testMateHitEquivalent();
    void testSingletonLessThanCigarOpLen();
    void testSingletonLessThanCigarOpCode();
    void testSingletonLessThanCigarLength();
    void testSingletonLessThanInsertID();
    void testSingletonEquivalant();
};

#endif