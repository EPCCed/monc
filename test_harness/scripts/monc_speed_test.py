#!/usr/bin/env python

# work out and output the speed of each simulation
# based on the contents of the monc stdout files of test_harness simulations
#   Uses bash script: monc_speed_check
#
# This is based on monc_kgo_bit_compare.py and should be incorporated at sometime
# For now, keep separate

import sys
import os
import subprocess as sub
import getopt

#import local modules
import dictionary_setup
import list_setup
import test_harness_setup

suite = ''
checkpoint_test = True
compiler = 'cray'
main_kgo_dir =  'dummy'
casim_kgo_dir = 'dummy'
new_monc_datadir = None

if len(sys.argv) > 1 :
    options, args = getopt.getopt(sys.argv[1:],"ht:c:p:",["test_suite=", "compiler=", "path="])

    for opt, arg in options:
        if opt == '-h' :
            print('Usage of monc_speed_test.py')
            print('monc_speed_test.py -t <standard, main_component, ecse, casim_socrates> -c <cray, gnu> ')
            print('-t or --test_suite is the suite type, i.e. component or standard')
            print('-c or --compiler is the compiler, i.e. cray or gnu, default is cray')
            print('-p or --path is the path to the test data (e.g., casedir/monc_stdout)')
            sys.exit()
        if opt in ("-t", "--test_suite"):
            suite = arg
        if opt in ('-c', '--compiler'):
            compiler = arg 
        if opt in ('-p', '--path'):
            new_monc_datadir = arg

else :
    print('No test-suite name provided')
    print('STOP - Must provide the suite, i.e. -t <standard, main_component, ecse, casim_socrates>' )
    print('Usage of monc_speed_test.py')
    print('monc_speed_test.py -t <standard, main_component, ecse, casim_socrates> -c <cray, gnu>')
    print('-t or --test_suite is the suite type, i.e. component or standard')
    print('-c or --compiler is the compiler, i.e. cray or gnu, default is cray')
    print('-p or --path is the path to the test data (e.g., casedir/monc_stdout)')
    sys.exit()

# read in lists and dictionaries, which are
# required to setup and run the test_harness
dict_list = test_harness_setup.dict_list_read(suite)
procs = dict_list[0]
test_case = dict_list[1]

if new_monc_datadir is None:
    new_monc_datadir = input("Please provide full path for directory with test data (monc_stdout): ")

for case in test_case :
    if  suite == 'main_component' : 
        print(case)
        test_list = list_setup.main_test_names(case)
        monc_dict = dictionary_setup.monc_main(case, main_kgo_dir, new_monc_datadir)
    elif suite == 'ecse':
        test_list = list_setup.ecse_test_names(case)
        monc_dict = dictionary_setup.monc_ecse(case, main_kgo_dir, new_monc_datadir)
    elif suite == 'casim_socrates' :
        test_list = list_setup.casim_test_names(case)
        monc_dict = dictionary_setup.monc_casim(case, casim_kgo_dir, new_monc_datadir)
    else:
        print('Unrecognised test_suite.  Re-try.')
        sys.exit()

    for test_idx, test_run in enumerate(test_list) :
        
        #runtime_calc.write(case, test_run, monc_dict)
        fn_36pe =  monc_dict['new_dir']+'output_'+test_run+monc_dict['no_procs'][0]+'*'
        fn_36_file = monc_dict['new_dir']+'times_'+test_run+monc_dict['no_procs'][0]+'.txt'
        p=sub.call(['monc_speed_check', fn_36pe, fn_36_file])
        
        text_file = open(fn_36_file, "r")
        lines = [x.strip('\n') for x in text_file.readlines()]
        lines_int = [int(x) for x in lines]
        
        time_36_pe = sum(lines_int)
        
        fn_72pe =  monc_dict['new_dir']+'output_'+test_run+monc_dict['no_procs'][1]+'*'
        fn_72_file = monc_dict['new_dir']+'times_'+test_run+monc_dict['no_procs'][1]+'.txt'
        p=sub.call(['monc_speed_check', fn_72pe, fn_72_file])

        text_file = open(fn_72_file, "r")
        lines = [x.strip('\n') for x in text_file.readlines()]
        lines_int = [int(x) for x in lines]
        
        time_72_pe = sum(lines_int)

        print(test_run, time_36_pe, time_72_pe)





