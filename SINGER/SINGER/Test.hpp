//
//  Test.hpp
//  SINGER
//
//  Created by Yun Deng on 3/17/23.
//

#ifndef Test_hpp
#define Test_hpp

#include <stdio.h>
#include "Sampler.hpp"
#include "Trace_pruner.hpp"
#include "Normalizer.hpp"

void test_read_arg();

void test_parsimony_pruner();

void test_pruner_efficiency();

void test_trace_pruner();

void test_iterative_start();

void test_fast_iterative_start();

void test_fast_larger_iterative_start();

void test_terminal_sampling();

void test_fast_terminal_sampling();

void test_internal_sampling();

void test_fast_internal_sampling();

void test_fast_larger_internal_sampling();

void test_optimal_ordering();

void test_load_vcf();

void test_no_recomb();

void test_no_mut();

void test_normalizer();

void test_read_muts();

void test_resume();

void test_fast_resume();

void test_resume_fast_larger_internal_sampling();

void test_tsp();

void test_debug_resume();

void test_african_dataset();

void test_start_african_sampling();

void test_resume_african_dataset();

#endif /* Test_hpp */
