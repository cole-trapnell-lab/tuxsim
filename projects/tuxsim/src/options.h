#ifndef OPTIONS_H
#define OPTIONS_H
/*
 *  options.h
 *  tuxsim
 *
 *  Created by Cole Trapnell on 4/14/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#define PACKAGE_VERSION "INTERNAL"
#endif

#include "common.h"

int parse_options(int argc, char** argv);
void print_usage();
#endif
