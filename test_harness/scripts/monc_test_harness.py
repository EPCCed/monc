#!/usr/bin/env python


# Modified to run only with Python3


# vn1.0 test-suite 
# This has been developed to automate the component testing of MONC, CASIM and SOCRATES.
# This is done by running a suite of tests that incrementally increase the complexity of 
# of the simulation by switching on components within the case. The component 
# testing is invoked by setting suite = 'component' and this will run the 
# configurations in the test_suite directory
# 
# As well as the component testing, this script can also run the standard cases
# so that a user can check that code change has not broken any of the standard cases.
# This is invoked by setting suite = 'standard'. This will use the configurations in 
# the testcases directory, but will copy the configurations and io configurations to
# the test_suite directory, so that the testcases remain stand alone.
#
# This script is an example of funtionality required for the development of a 
# Rose suite to automate the testing of the MONC. The script below is specific to MONSOON
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
#
# This test harness is a combination of the monc_main_test_suite and monc_casim_test_suite 
#
# Adhill - 071218

import subprocess as sub
import os
import sys
import getopt
from util import strtobool
from shutil import copyfile,copy,copytree,rmtree
import fileinput

# python version check
#   - however, this script will fail on a syntax error for the 1st f-string if using Python 2.
if sys.version_info.major != 3:
    print(' STOP -- The test_harness scripts now work with Python 3 only.')
    print('    On Monsoon, try: module load scitools')
    print('    On ARCHER2, try: module load cray-python')
    sys.exit()

# local scripts
import test_harness_setup # this function is in test_harness_setup.py
from util import yes_or_no, replaceInFile, replaceLineInFile, appendToFile

# Usage message
def print_usage_message():
    print('\n-------------------------------------------------------------------------')
    print('Usage of monc_test_harness.py\n')
    print('monc_test_harness.py -t <standard, main_component, ecse, casim_socrates, casim_aerosol_processing> -b <True or False> -m <default, cdt> -o <debug, safe, high> -c <cray, gnu> -s <MO-XC40, ARCHER2> -a <project-account>, -r <True or False>, -l <optional string label>, -i <?option1=value1?option2=value2?...>\n')
    print('Options:')
    print('    -t or --test_suite,    DEFAULT: main_component')
    print('        the suite type, i.e. standard, main_component, ecse, casim_socrates, or casim_aerosol_processing')
    print('        N.B.: The standard suite uses existing directories in the local branch, while')
    print('              main_component and casim_socrates use new directories under test_harness.')
    print('              These have not yet been configured to run under Slurm (ARCHER2).\n')
    print('    -b or --build,         DEFAULT: False')
    print('        whether to build MONC, True or False\n')
    print('    -r or --run,           DEFAULT: False')
    print('        whether to run the test harness, True or False\n')
    print('    -m or --modules,       DEFAULT: default')
    print('        the modules to load, i.e. default or cdt, where cdt is latest modules from')
    print('        Cray Development tools (Monsoon ONLY). Also, on ARCHER2: md1.0, md2.0, 2.0-cce20 \n')
    print('    -c or --compiler,      DEFAULT: cray')
    print('        the compiler, i.e. cray or gnu\n')
    print('    -o or --opt_level,     DEFAULT: high')
    print('        the level of optimisation and compiler flags, i.e. debug, safe, high\n')
    print('    -s or --system,        DEFAULT: MO-XC40')
    print('        the system to run the test suite on, MO-XC40 or ARCHER2\n')
    print('    -a or --account,       DEFAULT: ""')
    print('        REQUIRED for ARCHER2, the name of the project account to charge\n')
    print('    -l or --label,         DEFAULT: ""')
    print('        Completely optional label appended to test_harness run directory.  Useful for repeated tests with')
    print('        the same test_harness options, but with potentially different model options or simple retests.\n')
    print('    -i or --input_option,  DEFAULT: ""')
    print('        A list of option/value pairs to be added to the configuration files for all cases in this')
    print('        test_suite.  These are given as a single string, beginning with and delimited by the ? symbol.')
    print('        For example, we can ensure that threading is disabled in the IO server and that the flux_budget')
    print('        component is disabled in all cases by supplying:')
    print('            -i  ?l_thoff=.true.?flux_budget_enabled=.false.')
    print('        Any number of these can be supplied.  They are not checked for validity until used by the model.')
    print('        Quotation marks are not retained when passed in this way.\n')
    print('default settings =')
    print('    monc_test_harness.py -t main_component -b False -r False -m default -c cray -o high -s MO-XC40')
    print('-------------------------------------------------------------------------\n')
    sys.exit()


# default script settings
suite = 'main_component'
build = False
compiler = 'cray'
modules = 'default'
run = False
system = 'MO-XC40'
batch = 'PBS'
pes_per_node = 36 
opt_level = 'high'
account = ''
label = ''
jlabel =''    # This is NOT user input, just a variable that can be used to set a job name that needs 
              # to be include the `label` tag, if it is used.
input_option=[]


# Read in user-supplied options
if len(sys.argv) > 1 :
    options, args = getopt.getopt(sys.argv[1:],"ht:c:s:b:r:o:m:a:l:i:",
                                               ["test_suite=", "build=", "compiler=", 
                                                "system=", "run=", "opt_level=",
                                                "modules=", "account=", "label=",
                                                "input_option=",])

    for opt, arg in options:
        if opt == '-h' :
            print_usage_message()
        if opt in ("-t", "--test_suite"):
            suite = arg
        if opt in ("-s", "--system"):
            system = arg   
        if opt in ("-b", "--build"):
            build = strtobool(arg)
        if opt in ('-m', '--modules') :
            modules = arg
        if opt in ('-c', '--compiler'):
            compiler = arg
        if opt in ('-r', '--run'):
            run =  strtobool(arg)  
        if opt in ('-o', '--opt_level'):
            opt_level = arg
        if opt in ('-a', '--account'):
            account = arg
        if opt in ('-l', '--label'):
            label = arg
        if opt in ('-i', '--input_option'):
            input_option = arg.split("?")

else :
    print('No options supplied.  Running default options:')
    print('    monc_test_harness.py -t main_component -b False -r False -m default -c cray -o high -s MO-XC40\n')
    

# Begin test harness set-up

# Create short reference to modules for jobname
mod_mapping = {
    'default':     'd',
    'cdt':         'cdt',
    'md1.0':       'c15',
    'md2.0':       'c16',
    'md2.0-cce20': 'c20' }
mod_code = mod_mapping[modules]


# Verify this script is being run from a MONC main directory
if os.path.exists(os.path.join(os.getcwd(), 'test_harness')):
    pass
else:
    print('The test_harness directory does not exist in this location.')
    print('STOP - Try running from a MONC main directory.')
    sys.exit()

# System-specific parameters
if system == 'MO-XC40' :
    pes_per_node = 36
    batch = 'PBS'
    command = 'qsub'
    qsub_script_template = 'test_harness/submonc_template.pbs'
elif any([system == s for s in ["ARCHER2", "ecse"]]) :
    pes_per_node = 128
    batch = 'SLURM'
    command = 'sbatch'
    qsub_script_template = 'test_harness/submonc_template.sb'
    if account == '':
        print('Attempting to configure for ARCHER2.')
        print('STOP - Must provide project account name, e.g. -a n02-REVCON')
        sys.exit()
    if suite == 'standard':
        print('Attempting to test the "standard" cases on ARCHER2.')
        print('STOP - These are not yet configured for ARCHER2.')
        sys.exit()
else :
    pes_per_node = input("Please enter the number of cores per node: ")


# default io config directory, use the standard monc io configs
#    NOTE: might wish to isolate this (and included files) under the test_harness directory,
#          as users might have modified these (very common)
ioserver_cfg_dir = 'io/io_cfg_files/'
ioserver_cfg_file=ioserver_cfg_dir+'data_write_1file.xml' 
    
# read in lists and dictionaries, which are
# required to setup and run the test_harness
dict_list = test_harness_setup.dict_list_read(suite)
procs = dict_list[0]
testcases = dict_list[1]
test_mcf = dict_list[2]

# Create directory for this specific run, if needed
master_dir = system + '_' + suite + '_' + compiler + '_' + modules + '_' + opt_level
if len(label) > 0:
    master_dir = master_dir + '_' + label
    jlabel = label + '_'
compiler_dir = 'test_harness/' + master_dir
if not os.path.exists(compiler_dir):
    os.makedirs(compiler_dir)

# New directory for build logs and executable
build_dir = compiler_dir + '/build_dir/'
if not os.path.exists(build_dir):
    os.makedirs(build_dir)

# Build MONC, if requested. 
if build :
    print(f'\n  [WARNING] Building MONC executable.  Removing existing.')
    try:
        os.remove('./build/bin/monc_driver.exe')
    except:
        print(f'\n  [WARNING] No existing executable...  Continuing...)')
    p=sub.call(['monc_build_cray.sh', compiler, system, opt_level, modules, build_dir])
else:
    print(f'\n  [WARNING] This is not a fresh build.  Copying existing executable to {build_dir}.\n')
    if not yes_or_no(f'Would you like to continue with the existing executable?'):
        print('Exiting')
        sys.exit()


# In all cases, copy the executable to the build_dir if it exists   
if os.path.isfile('./build/bin/monc_driver.exe'):
    # Orignal executable location
    build_exec = './build/bin/monc_driver.exe'
    
    # Copy the executable to the new location
    copy(build_exec,build_dir+'monc_driver.exe')

    # Copy the build logs
    if os.path.exists(build_dir + '.fcm-make/'):
        rmtree(build_dir + '.fcm-make/')
    if os.path.exists(build_dir + 'preprocess/'):
       rmtree(build_dir + 'preprocess/')
    copytree('.fcm-make/', build_dir + '.fcm-make/')    # copytree was slow in testing and requires
    copytree('preprocess/', build_dir + 'preprocess/')  # rmtree to overwrite existing

    # Copy the version file if it exits (fcm output and the module info)
    if os.path.isfile('./version'):
        copy('./version',build_dir+'version')
else:
    print('\nNo executable found.\n')
    if not yes_or_no(f'Would you like to continue creating directories and files anyway?'):
        print('Exiting')
        sys.exit()

# Create an empty directory to later store bit comparison information
bc_dir = compiler_dir + '/bit_comparison/'
if not os.path.exists(bc_dir):
    os.makedirs(bc_dir)

# Set up the PBS out directory (same name on ARCHER2)
# Contains the job -oe files - reduces clutter in main directory.
# Since this redirection throws off the naming, a record of the 
# $PBS_JOBIDs are printed to the monc_stdout.
# Alternatively, we could change the submission directory, but
# since the paths are relative everywhere, we'd have to change 
# many things...
pbs_dir = compiler_dir + '/pbs_dir/'
if not os.path.exists(pbs_dir):
    os.makedirs(pbs_dir)

## This script makes the following directories, if they do not exist
##
# this directory contains the configs that are run
# the files are built using the base configs in bubble, stratus
# directories. These will all be the same except the xml file
tmp_config_dir =compiler_dir+'/tmp_config_dir/'
#tmp_ioserver_cfg_dir = compiler_dir+'/tmp_ioserver_cfg_dir/'
# Submission scripts are saved in this directory
tmp_qsub_dir = compiler_dir+'/tmp_qsub_dir/'
# Netcdf output is saved in this directory
nc_output_dir = compiler_dir+'/ncfiles/'
# Checkpoint output is save in this directory
checkpoint_output_dir = compiler_dir+'/chk_point_dump/'
# the standard output is stored here. This is required for 
# debug and can be used for control of the test harness
out_dir = compiler_dir+'/monc_stdout/'

if not os.path.exists(tmp_config_dir):
    os.makedirs(tmp_config_dir)
    
#if not os.path.exists(tmp_ioserver_cfg_dir):
#    os.makedirs(tmp_ioserver_cfg_dir)

if not os.path.exists(tmp_qsub_dir):
    os.makedirs(tmp_qsub_dir)
    
if not os.path.exists(nc_output_dir):
    os.makedirs(nc_output_dir)
    
if not os.path.exists(checkpoint_output_dir):
    os.makedirs(checkpoint_output_dir)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
for tc in testcases :
    # Setup base dir name for the suite testcases 
    if suite == 'main_component' :
        base_cfg_dir = 'test_harness/monc_main/'+tc+'/'
    elif suite == 'casim_socrates' :
        base_cfg_dir = 'test_harness/monc_casim_socrates/'+tc+'/'
    elif suite == 'ecse':
        base_cfg_dir = 'test_harness/monc_ecse/'+tc+'/'
    elif suite == 'casim_aerosol_processing' :
        base_cfg_dir = 'test_harness/casim_aerosol_processing/'+tc+'/'
    elif suite == 'standard' :
        base_cfg_dir = 'testcases/'+tc+'/'
    else :
        sys.exit('suite name not recognised - STOP')
    
    for cfg in test_mcf[tc]:
        print('')

        for pe in procs : 

      ### The standard case doesn't actually do pe tests, so it shouldn't be in this loop or it could remain here with a conditional break after the submission section.
            if suite == "standard" :
                # 1. set up config file and write to the temporary location (do not really
                #    need to do this, so might be worth removing the copy to tmp)
                config_file = base_cfg_dir+cfg+'.mcf'
                new_cfg_file = tmp_config_dir+cfg+'.mcf'
                copyfile(config_file, new_cfg_file)
                
                # 2 use default script from the testcase
                qsub_scriptfile_name = base_cfg_dir+'/submonc_scripts/'+cfg+'.pbs'
                   # -- NEEDS an ARCHER2 ALTERNATIVE, PERHAPS A CONVERSION ROUTINE --
                   # 1. Most of the script can be copied.
                   # 2. Need to replace directives.
                   #    - might need to read existing select/walltimes to get appropriate combinations.
                   # 3. Need to replace the WORKDIR items (as in test_harness template).
                   # 4. There is NO ck_progress handling of the standard case.
                new_scriptfile = qsub_scriptfile_name
            else :
                # 1. set up config file and write to the temporary location
                config_file = base_cfg_dir+cfg
                new_cfg_file = tmp_config_dir+cfg+'_'+str(pe)
                # 2. set up the netcdf output file name and add to the ioserver_cfg file
                output_data = nc_output_dir+cfg+'_'+str(pe)
                # 3. setup the checkpoint file name in the config file
                checkpoint_file=checkpoint_output_dir+cfg+'_'+str(pe)+'.nc'
                #replaceInFile(new_cfg_file, 'add_checkpoint_file', checkpoint_file, new_cfg_file)
                replaceInFile(config_file, 'add_checkpoint_file', checkpoint_file, new_cfg_file)
                # 4. set up diagnostic file name in the .mcf (all diags written to one file)
                replaceInFile(new_cfg_file, 'add_diagnostic_file', output_data+'_diagnostics.nc', new_cfg_file)                               
                # 5. set up qsub submission file in the test-harness tmp directory. 
                #    There is a submission script for each run
                new_scriptfile = tmp_qsub_dir+'submonc_'+cfg+'_'+str(pe) 
                number_of_nodes = pe // pes_per_node
                mod_pes = pe % pes_per_node
                if mod_pes > 0 :
                    number_of_nodes = number_of_nodes + 1
                print('Job pes: {!r}, pes_per_node: {!r}, mod_pes: {!r}, number_of_nodes: {!r}'.format(pe, pes_per_node, mod_pes, number_of_nodes))
                if batch == 'PBS':
                    replaceInFile(qsub_script_template, 'select=1', 'select='+str(number_of_nodes), new_scriptfile)
                elif batch == 'SLURM':
                    replaceInFile(qsub_script_template, '--nodes=1', '--nodes='+str(number_of_nodes), new_scriptfile)

                replaceInFile(new_scriptfile, 'add_np', str(pe), new_scriptfile)
                replaceInFile(new_scriptfile, 'add_scriptname', new_scriptfile, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_testcase_mcf', new_cfg_file, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_jobname', jlabel + suite[0] +'_'+ compiler[0] +'_'+ mod_code + '_' + opt_level[0] +'_'+ cfg +'_'+ str(pe), new_scriptfile)
                replaceInFile(new_scriptfile, 'add_stdout_dirname', out_dir, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_checkpoint_dirname', checkpoint_output_dir, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_testsuite_jobname', cfg+'_'+str(pe)+'_' , new_scriptfile)
                replaceInFile(new_scriptfile, 'add_build_exec','./'+build_dir+'monc_driver.exe', new_scriptfile)
                replaceInFile(new_scriptfile, 'add_pbs_dir', pbs_dir, new_scriptfile)
                if system == 'ARCHER2':
                    replaceInFile(new_scriptfile, 'add_account', account, new_scriptfile)
                    replaceInFile(new_scriptfile, 'add_modules', 'archer_'+compiler+'_mod', new_scriptfile)
                # 6. Append user supplied options_database entries to all configs
                if len(input_option) > 0:
                    appendToFile(new_cfg_file, '\n\n# Additional options added by monc_test_harness.py\n')
                    for ode in input_option:
                        appendToFile(new_cfg_file, f'{ode}\n')

### This might not be needed because we've already set up to conform to the 36/72 core tests.
### However, that means that we aren't using mutliple nodes, which might be part of the intention of the test...
#                if system == 'ARCHER2':
#                    replaceLineInFile(new_scriptfile, 'moncs_per_io_server', 'moncs_per_io_server=15')


            # Finally determine whether to run the test harness or just set up the directories and files        
            if run :


# Might want to replace this with a call to ck_progress.sh (-a -o -c "specific_directory_name")
# I haven't tested this yet.  Not sure that I want to do it this way. -TRJ
#    p=sub.call(['ck_progress.sh', '-a', '-o', '-c', master_dir])

                print('submitting '+new_scriptfile+' to HPC')
                p=sub.call([command, new_scriptfile ])

            else :
               print('mcf '+new_cfg_file+' and script '+new_scriptfile+' setup but not run')
               print(command+' '+new_scriptfile)




        
        
        
        

