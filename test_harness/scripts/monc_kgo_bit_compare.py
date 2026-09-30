#!/usr/bin/env python


# This script has been updated to work only with Python 3.

# Beta test-suite (see ticket #65) 
# This has been developed to automate the component testing of MONC
# The attached script is an example of funtionality required for the development of a 
# Rose suite to automate the testing of the MONC. The script below is specific to MONSOON
# The script below is used to analyse the output from monc_test_harness.py and compare with 
# the LEM, where appropriate.
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
# Adhill - 040516

import sys
import os
import getopt

# python version check
if sys.version_info.major != 3:
    print(' STOP -- The test_harness scripts now work with Python 3 only.')
    print('    On Monsoon, try: module load scitools')
    print('    On ARCHER2, try: module load cray-python')
    sys.exit()

#import local module
import bit_compare
import dictionary_setup
import list_setup
import test_harness_setup
from util import yes_or_no

# Usage message
def print_usage_message():
    print('\n-------------------------------------------------------------------------')
    print('Usage of monc_kgo_bit_compare.py')
    print('monc_kgo_bit_compare.py -t <standard, main_component, ecse, casim_socrates, casim_aerosol_processing> -c <cray, gnu> ')
    print('Options:')
    print('    -t or --test_suite,  DEFAULT: ""')
    print('        the suite type, i.e. standard, main_component, ecse, casim_socrates, casim_aerosol_processing\n')
    print('    -c or --compiler     DEFAULT: ""')
    print('        the compiler, i.e. cray or gnu\n')
    print('    -o or --opt_level    DEFAULT: ""')
    print('        optimisation level, i.e. high, safe, debug\n')
    print('    -k or --kgo_path,    DEFAULT: undefined')
    print('        absolute or relative path to KGO data directory, e.g., ../ncfiles')
    print('        if not provided, user will be prompted to enter this information\n')
    print('    -n or --new_path,    DEFAULT: undefined')
    print('        absolute or relative path to New data directory, e.g., ../ncfiles')
    print('        if not provided, user will be prompted to enter this information\n')
    print('NOTE: the compiler and opt_level options are only printed to the output file. They are not checked.')
    print('\nWhen running, the program will request user to input the locations of the two datasets to be compared.')
    print('\nTypically, we like to run the KGO test from the test_harness/<master_dir>/bit_comparison directory.')
    print('-------------------------------------------------------------------------\n')
    sys.exit()

# default script settings
suite = ''
opt_level = ''
compiler = ''
input_kgo = None
new_monc_datadir = None
plot_all = False

# Read in user-supplied options
if len(sys.argv) > 1 :
    options, args = getopt.getopt(sys.argv[1:],"ht:c:o:k:n:a",
                                 ["test_suite=", "compiler=", "opt_level=",
                                  "kgo_path=", "new_path="])

    for opt, arg in options:
        if opt == '-h' :
            print_usage_message()
        if opt in ('-t', '--test_suite'):
            suite = arg
        if opt in ('-c', '--compiler'):
            compiler = arg 
        if opt in ('-o', '--opt_level'):
            opt_level = arg 
        if opt in ('-k', '--kgo_path'):
            input_kgo = arg
        if opt in ('-n', '--new_path'):
            new_monc_datadir = arg
        if opt == '-a':
            plot_all = True

else :
    print( 'No test-suite name provided')
    print( 'STOP - Must provide the suite, i.e. -t <standard or main_component or casim_socrates>' )
    print_usage_message()


# Check whether this script is being run from the bit_comparison directory
if os.getcwd().split('/')[-1] == 'bit_comparison':
    pass
else:
    print('Typically, we like to run the KGO test from the test_harness/<master_dir>/bit_comparison directory.')
    if not yes_or_no(f'Would you like to run from {os.getcwd()} instead?'):
        print('Exiting')
        sys.exit()


# read in lists and dictionaries, which are
# required to setup and run the test_harness
dict_list = test_harness_setup.dict_list_read(suite)
procs = dict_list[0]
test_case = dict_list[1]


# Request kgo directory, if not already provided
if input_kgo is None:
    input_kgo = input("Please provide full path to KGO directory: ")
if not os.path.isdir(input_kgo) :
    sys.exit(input_kgo+' does not exist, please check and try again - Exiting')


# Request new test data directory, if not already provided
if new_monc_datadir is None:
    new_monc_datadir = input("Please provide full path to directory with test data: ")
if not os.path.isdir(new_monc_datadir) :
    sys.exit(new_monc_datadir+' does not exist, please check and try again - Exiting')


# Convert to absolute paths, in case someone passes a '../ncfiles'
# Also ensures there's a trailing slash on the path.
input_kgo = os.path.join(os.path.abspath(input_kgo),'')
new_monc_datadir = os.path.join(os.path.abspath(new_monc_datadir),'')


# Open results file
monc_kgo_datadir = input_kgo
with open('bit_compare_results.txt', "w") as f:
    f.write('****************Bit comparison results from  '+suite+'  component testing**************'+"\n")
    f.write('opt_level: ' + opt_level + "\n")
    f.write('compiler:  ' + compiler + "\n")
    f.write('Known good MONC output = '+monc_kgo_datadir+"\n")
    f.write('New MONC output = '+new_monc_datadir+"\n")

# Work through the test cases
for case in test_case :
    with open('bit_compare_results.txt', "a") as f:
        f.write('********************************************************************************************'+"\n")
        f.write('******************************************'+case+'******************************************'+"\n")

    print(f'\n{case}')

    if  suite == 'main_component' : 
        test_list = list_setup.main_test_names(case)
        figure_title_list = list_setup.main_figure_names(case)
        monc_dict = dictionary_setup.monc_main(case, monc_kgo_datadir, new_monc_datadir)
    elif suite == 'ecse':
        test_list = list_setup.ecse_test_names(case)
        figure_title_list = list_setup.ecse_figure_names(case)
        monc_dict = dictionary_setup.monc_ecse(case, monc_kgo_datadir, new_monc_datadir)
    elif suite == 'casim_socrates' :
        test_list = list_setup.casim_test_names(case)
        figure_title_list = list_setup.casim_figure_names(case)
        monc_dict = dictionary_setup.monc_casim(case, monc_kgo_datadir, new_monc_datadir)
    elif suite == 'casim_aerosol_processing' :
        test_list = list_setup.casim_aeroproc_test_names(case)
        figure_title_list = list_setup.casim_aeroproc_figure_names(case)
        monc_dict = dictionary_setup.monc_casim_aeroproc(case, monc_kgo_datadir, new_monc_datadir)

    # do bit comparison check on MONC and maybe LEM
    for test_idx, test_run in enumerate(test_list) :
        print(test_run)
        with open('bit_compare_results.txt', "a") as f:
            f.write('=========================================================================================='+"\n")
        bit_compare_success = bit_compare.monc_calc(case, test_run, monc_dict)
        
        if not bit_compare_success[0] or plot_all:  
            print('      !! plotting bit comparison from monc to assess issue !!')
            bit_compare.monc_plot('bit_compare', case, test_run, figure_title_list[test_idx], monc_dict) 

            if bit_compare_success[1] :
                bit_compare.monc_plot('kgo_compare_36pe', case, test_run, figure_title_list[test_idx], monc_dict)  
                bit_compare.monc_plot('kgo_compare_72pe', case, test_run, figure_title_list[test_idx], monc_dict)
                


