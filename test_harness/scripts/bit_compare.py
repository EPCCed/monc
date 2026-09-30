#!/usr/bin/env python

# Beta test-suite (see ticket #65)
# Code for calculating and plotting the bit comparison tests 
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
# Adhill - 040516

import netCDF4
import numpy as np
import matplotlib.pyplot as p
import os
import sys
import subprocess as sub
#from prettytable import PrettyTable

#import local modules
import get_data
import list_setup

# Const
big_err = 999999999

def monc_calc(testcase, test_run, model_dict) :
    
    Bit_compare_threshold = 0.0

    # list of 3d fields to compare
    # list of 3d fields to compare

# Remind me, why does rce exclude moist variables? 
    if testcase == 'stratus' or testcase == 'shallow_convection' :
        fields_3d = ['u', 'v', 'w', 'q_vapour', 'q_cloud_liquid_mass']
    else :
        fields_3d = ['u', 'v', 'w' ]


    # Check whether decomposition failure is expected; create related messaging
    #   (cases that include a global mean or sum, e.g. from damping or forcing)
    decomp_comm = '[decomp failure expected]' if test_run in list_setup.expected_decompostion_failures else '[decomp success expected]'

    # Annotate file with case under consideration
    with open('bit_compare_results.txt', "a") as f:
        f.write(f'  Case:  {test_run}  {decomp_comm}\n')

    # Work out the known good output (kgo) difference
    #
    # construct filename using model_dict information
    # KGO first
    new_file_test = True
    kgo_file_test = True
    kgo_fn_36pe =  model_dict['kgo_dir']+test_run+model_dict['no_procs'][0]+model_dict['files'][0]
#    print(kgo_fn_36pe)
    if (not os.path.exists(kgo_fn_36pe)) :
        kgo_file_test = False
        print("Missing: ",kgo_fn_36pe )
    kgo_fn_72pe = model_dict['kgo_dir']+test_run+model_dict['no_procs'][1]+model_dict['files'][0]
#    print(kgo_fn_72pe)
    if (not os.path.exists(kgo_fn_72pe)) :
        kgo_file_test = False
        print("Missing: ", kgo_fn_72pe)
    # New ouput 
    new_fn_36pe =  model_dict['new_dir']+test_run+model_dict['no_procs'][0]+model_dict['files'][0]
#    print(new_fn_36pe)
    if (not os.path.exists(new_fn_36pe)) :
        new_file_test = False
        print("Missing: ", new_fn_36pe)
    new_fn_72pe = model_dict['new_dir']+test_run+model_dict['no_procs'][1]+model_dict['files'][0]
#    print(new_fn_72pe)
    if (not os.path.exists(new_fn_72pe)) :
        new_file_test = False
        print("Missing: ", new_fn_72pe)

    # check files exist
    if not new_file_test :
        with open('bit_compare_results.txt', "a") as f:
            f.write(test_run+' bit compare **FAILED** - file missing, continue to next test'+"\n")
        return True, kgo_file_test  
    else : # new files exist, so carry on with bit comparison
        # open both KGO files (if they exist)
        # script will still do bit compare test on new files even if
        # the kgo does not exist
        if kgo_file_test :
            kgo_36pe_open=netCDF4.Dataset(kgo_fn_36pe,'r') 
            kgo_72pe_open=netCDF4.Dataset(kgo_fn_72pe,'r')
 
        # open both NEW test files
        new_36pe_open=netCDF4.Dataset(new_fn_36pe,'r') 
        new_72pe_open=netCDF4.Dataset(new_fn_72pe,'r')
     
        # Work out the difference the new output on different decompositions 
        all_kgo_decomp_success = True
        all_new_decomp_success = True
        all_kgo_new_compare_success = True

        # Loop over relevant fields to compare
        for field in fields_3d :

            # Annotate field to compare
            with open('bit_compare_results.txt', "a") as f:
                f.write('Field = '+field+'\n')

            kgo_decomp_success = True
            new_decomp_success = True
            kgo_new_compare_success = True

            if kgo_file_test :
                # Extract field from the KGO file
                kgo_w_36pe = kgo_36pe_open.variables[field][:,:,:,:]
                kgo_w_72pe = kgo_72pe_open.variables[field][:,:,:,:]

            # Extract field from the New output
            new_w_36pe = new_36pe_open.variables[field][:,:,:,:]
            new_w_72pe = new_72pe_open.variables[field][:,:,:,:]
        
            # work out the difference between the processor decompositions
            # KGO
            if kgo_file_test :
                try:
                    kgo_diff_pe = kgo_w_36pe - kgo_w_72pe
                # This will catch any error in the calculation, but this is often due to a mismatch in array shape
                except Exception as err:
                    print(f"Unexpected err={err}, type(err)={type(err)}")
                    print(f"kgo_w_36pe.shape={kgo_w_36pe.shape}, kgo_w_72pe.shape={kgo_w_72pe.shape}")
                    with open('bit_compare_results.txt', "a") as f:
                            f.write(f"Unexpected err={err}, type(err)={type(err)}\n")
                            f.write(f"kgo_w_36pe.shape={kgo_w_36pe.shape}, kgo_w_72pe.shape={kgo_w_72pe.shape}\n")
                    kgo_diff_pe = big_err

                if (np.max(abs(kgo_diff_pe)) > Bit_compare_threshold) :
                    decomp_success = False 
                    kgo_decomp_success = False
            # NEW 
            try:
                new_diff_pe = new_w_36pe - new_w_72pe
            # This will catch any error in the calculation, but this is often due to a mismatch in array shape
            except Exception as err:
                print(f"Unexpected err={err}, type(err)={type(err)}")
                print(f"new_w_36pe.shape={new_w_36pe.shape}, new_w_72pe.shape={new_w_72pe.shape}")
                with open('bit_compare_results.txt', "a") as f:
                    f.write(f"Unexpected err={err}, type(err)={type(err)}\n")
                    f.write(f"new_w_36pe.shape={new_w_36pe.shape}, new_w_72pe.shape={new_w_72pe.shape}\n")
                new_diff_pe = big_err


            # Annotate file with test under consideration
            with open('bit_compare_results.txt', "a") as f:
                f.write(f'Decomposition comparison:\n')

            # Report on decomposition success            
            if (np.max(abs(new_diff_pe)) > Bit_compare_threshold) :  # decomp failure (difference) in new MONC output
                new_decomp_success = False
                if kgo_file_test :          
                    if kgo_decomp_success :  # decomp success (match) in KGO
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' **FAILED** (but KGO Passed) -ATTENTION: new decomposition failure-'+"\n")
                            f.write('Field = '+field+' -   Max new decomp diff ='+str(np.max(abs(new_diff_pe)))+"\n")
                    else :                   # decomp failure (difference) in KGO
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' **FAILED** (KGO also Failed)'+"\n")
                            f.write('Field = '+field+' -   Max KGO decomp diff ='+str(np.max(abs(kgo_diff_pe)))+"\n")
                            f.write('Field = '+field+' -   Max new decomp diff ='+str(np.max(abs(new_diff_pe)))+"\n")
                else :
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' **FAILED** (NO KGO files to test)'+"\n")
                        f.write('Field = '+field+' -   Max new decomp diff ='+str(np.max(abs(new_diff_pe)))+"\n")
            else :                                                   # no decomp difference in new
                if kgo_file_test : 
                    if kgo_decomp_success :
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS (KGO also Passed)'+"\n")
                    else :
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS (but KGO Failed) -ATTENTION: new decomposition success; unusual...'+"\n")
                else :
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS (NO KGO files to test)'+"\n")


            # Annotate file with test under consideration
            with open('bit_compare_results.txt', "a") as f:
                f.write(f'New vs KGO comparison:\n')

            # Compare kgo and new (if kgo exists)
            if kgo_file_test :
                try:
                    diff_kgo_new_36pe = new_w_36pe - kgo_w_36pe

                # This will catch any error in the calculation, but this is often due to a mismatch in array shape
                except Exception as err:
                    print(f"Unexpected err={err}, type(err)={type(err)}")
                    print(f"new_w_36pe.shape={new_w_36pe.shape}, kgo_w_36pe.shape={kgo_w_36pe.shape}")
                    with open('bit_compare_results.txt', "a") as f:
                            f.write(f"Unexpected {err}, {type(err)}\n")
                            f.write(f"new_w_36pe.shape={new_w_36pe.shape}, kgo_w_36pe.shape={kgo_w_36pe.shape}\n")
                    diff_kgo_new_36pe = big_err

                try:
                    diff_kgo_new_72pe = new_w_72pe - kgo_w_72pe

                # This will catch any error in the calculation, but this is often due to a mismatch in array shape
                except Exception as err:
                    print(f"Unexpected err={err}, type(err)={type(err)}")
                    print(f"new_w_72pe.shape={new_w_72pe.shape}, kgo_w_72pe.shape={kgo_w_72pe.shape}")
                    with open('bit_compare_results.txt', "a") as f:
                            f.write(f"Unexpected {err}, {type(err)}\n")
                            f.write(f"new_w_72pe.shape={new_w_72pe.shape}, kgo_w_72pe.shape={kgo_w_72pe.shape}\n")
                    diff_kgo_new_72pe = big_err



                # Report difference to stdout
                kgo_new_compare_success = True
                if (np.max(abs(diff_kgo_new_36pe)) > Bit_compare_threshold) :
                    kgo_new_compare_success = False
                    print( "  KGO/New comparison difference 36pe", kgo_new_compare_success, field)
                if (np.max(abs(diff_kgo_new_72pe)) > Bit_compare_threshold ) :
                    kgo_new_compare_success = False
                    print( "  KGO/New comparison difference 72pe", kgo_new_compare_success, field)

                # Report difference to results file
                if (np.max(abs(diff_kgo_new_36pe)) > Bit_compare_threshold ) :  # compare 36s
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs kgo '+model_dict['no_procs'][0]+'pe '+test_run+' **FAILED** -ATTENTION: Science change!-'+"\n")  
                        f.write('Field = '+field+' -   Max New/kgo diff ='+str(np.max(abs(diff_kgo_new_36pe)))+"\n")
                else :                                                          # compare 72s if 36s match
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field = '+field+' - New '+model_dict['no_procs'][0]+'pe vs kgo '+model_dict['no_procs'][0]+'pe '+test_run+' SUCCESS'+"\n")

                if (np.max(abs(diff_kgo_new_72pe)) > Bit_compare_threshold ) :
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field = '+field+' - New '+model_dict['no_procs'][1]+'pe vs kgo '+model_dict['no_procs'][1]+'pe '+test_run+' **FAILED** -ATTENTION: Science change!-'+"\n") 
                        f.write('Field = '+field+' -   Max New/kgo diff ='+str(np.max(abs(diff_kgo_new_72pe)))+"\n")
                else :
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field = '+field+' - New '+model_dict['no_procs'][1]+'pe vs kgo '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS'+"\n")  
                with open('bit_compare_results.txt', "a") as f:
                    f.write('..........................................................................................'+"\n")
            # No kgo present
            else:
                with open('bit_compare_results.txt', "a") as f:
                    f.write('Field = '+field+' - New '+model_dict['no_procs'][1]+'pe vs kgo '+model_dict['no_procs'][1]+'pe '+test_run+' NO KGO DATA'+"\n") 
                    f.write('..........................................................................................'+"\n")

            # Update overall flag
            all_kgo_decomp_success = all_kgo_decomp_success and kgo_decomp_success
            all_new_decomp_success = all_new_decomp_success and new_decomp_success
            all_kgo_new_compare_success = all_kgo_new_compare_success and kgo_new_compare_success

        if kgo_file_test:            
            kgo_36pe_open.close()
            kgo_72pe_open.close()

        new_36pe_open.close()
        new_72pe_open.close()    

        print( f'    Note: {decomp_comm}')
        print( f'    Result: all_kgo_decomp_success: {all_kgo_decomp_success}, all_new_decomp_success: {all_new_decomp_success}, all_kgo_new_compare_success: {all_kgo_new_compare_success}')

        all_success = all([all_kgo_decomp_success,all_new_decomp_success,all_kgo_new_compare_success])
        return all_success, kgo_file_test  
      
def lem_calc(testcase, test_run, model_dict) :

    Bit_compare_threshold = 0.0

    for idx, file_name in enumerate(model_dict['files']) :
        # construct filename using model_dict information
        model_fn_36pe =  model_dict['kgo_dir']+test_run+model_dict['no_procs'][0]+file_name+'.nc'
        model_fn_72pe = model_dict['kgo_dir']+test_run+model_dict['no_procs'][1]+file_name+'.nc'
            
        # open both files
        f_36pe_open=netCDF4.Dataset(model_fn_36pe,'r') 
        f_72pe_open=netCDF4.Dataset(model_fn_72pe,'r')

        w_36pe = f_36pe_open.variables['W'][:,:,:]
        w_72pe = f_72pe_open.variables['W'][:,:,:]
            
        diff_pe = w_36pe - w_72pe

        if (np.max(abs(diff_pe)) > Bit_compare_threshold) :
            print( 'LEM '+test_run+' bit compare test also **FAILED**')
            return False
        else :
            if idx == len(model_dict['files'])-1 :
                print( 'LEM '+test_run+' bit compare test SUCCESS')
                return True

def monc_plot(testing, testcase, test_run, fig_title, monc_dict):
    
    #
    # plot the bit comparison results first
    #
    if testing == 'bit_compare':
        model_36Pe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][0]
        model_MultiPe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][1]
        label_prefix_36Pe = 'MONC 36 pe '
        label_prefix_72Pe = 'MONC 72 pe '
    elif testing == 'kgo_compare_36pe':
        model_36Pe_fn = monc_dict['kgo_dir']+test_run+monc_dict['no_procs'][0]
        model_MultiPe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][0]
        label_prefix_36Pe = 'KGO 36 pe '
        label_prefix_72Pe = 'NEW 36 pe '
    elif testing == 'kgo_compare_72pe': 
        model_36Pe_fn = monc_dict['kgo_dir']+test_run+monc_dict['no_procs'][1]
        label_prefix_36Pe = 'KGO 72 pe '
        label_prefix_72Pe = 'NEW 72 pe '
        model_MultiPe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][1]
        

    # this dictionary is used to output data from get_profile
    #it is overwritten on each loop
    if testcase == 'rce':
        monc_36Pe_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                          'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                          'vapour_mmr_mean': [], 'liquid_mmr_mean': [], 'rain_mmr_mean':[], 
                          'ice_mmr_mean': [], 'snow_mmr_mean':[],  'graupel_mmr_mean':[]  }
        get_data.profile_diagnostics(model_36Pe_fn, monc_dict['files'], monc_36Pe_data)
        
        # this dictionary is used to output data from get_profile 
        #it is overwritten on each loop
        monc_MultiPe_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                              'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                              'vapour_mmr_mean': [], 'liquid_mmr_mean': [],  'rain_mmr_mean':[], 
                              'ice_mmr_mean': [], 'snow_mmr_mean':[],  'graupel_mmr_mean':[] }
        get_data.profile_diagnostics(model_MultiPe_fn, monc_dict['files'], monc_MultiPe_data)
        
    else :    

        monc_36Pe_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                          'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                          'vapour_mmr_mean': [], 'liquid_mmr_mean': []}
        get_data.profile_diagnostics(model_36Pe_fn, monc_dict['files'], monc_36Pe_data)
        
        # this dictionary is used to output data from get_profile 
        #it is overwritten on each loop
        monc_MultiPe_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                              'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                              'vapour_mmr_mean': [], 'liquid_mmr_mean': []}
        get_data.profile_diagnostics(model_MultiPe_fn, monc_dict['files'], monc_MultiPe_data)

    z_plot =  monc_MultiPe_data['zn'][0][:]

    # Set up the velocity figure for the case 
    #
    velocity_fig = p.figure(figsize=(15,10))
    velocity_fig.suptitle(fig_title, fontsize=16 )
    ax1 = velocity_fig.add_subplot(231)
    ax2 = velocity_fig.add_subplot(232)
    ax3 = velocity_fig.add_subplot(233)
    ax4 = velocity_fig.add_subplot(234)
    ax5 = velocity_fig.add_subplot(235)
    ax6 = velocity_fig.add_subplot(236)

    for file_index, file_numbers in enumerate(monc_dict['files']) :
        ax1.plot(monc_36Pe_data['u_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index],label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
        ax1.plot(monc_MultiPe_data['u_mean'][0][file_index,:], z_plot[:]
                 , '--',color=monc_dict['time_color'][file_index],label=label_prefix_72Pe+monc_dict['t_labels'][file_index] )
            
        ax2.plot(monc_36Pe_data['v_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
        ax2.plot(monc_MultiPe_data['v_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            
        ax3.plot(monc_36Pe_data['w_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
        ax3.plot(monc_MultiPe_data['w_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
        
        ax4.plot(monc_36Pe_data['uu_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
        ax4.plot(monc_MultiPe_data['uu_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
        
        ax5.plot(monc_36Pe_data['vv_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
        ax5.plot(monc_MultiPe_data['vv_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
        
        ax6.plot(monc_36Pe_data['ww_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
        ax6.plot(monc_MultiPe_data['ww_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
    
    ax1.set_title('Mean U')
    ax2.set_title('Mean V')
    ax3.set_title('Mean W')
    ax4.set_title('Mean UU')
    ax5.set_title('Mean VV')
    ax6.set_title('Mean WW')
    ax1.set_ylabel('height (m)')
    ax4.set_ylabel('height (m)')
    ax1.set_xlabel('m s$^{-1}$')
    ax2.set_xlabel('m s$^{-1}$')
    ax3.set_xlabel('m s$^{-1}$')
    ax4.set_xlabel('m$^{2}$ s$^{-2}$')
    ax5.set_xlabel('m$^{2}$ s$^{-2}$')
    ax6.set_xlabel('m$^{2}$ s$^{-2}$')
        
    lgd = ax3.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
    
    p.savefig(testing+'_'+test_run+'velocity_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    p.close(velocity_fig)

    if testcase == 'bubble' :
        theta_fig = p.figure(figsize=(8,5))
        #theta_fig.suptitle(fig_title, fontsize=16 )
        ax1 = theta_fig.add_subplot(111)
        for file_index, file_numbers in enumerate(monc_dict['files']) :
            ax1.plot(monc_36Pe_data['theta_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax1.plot(monc_MultiPe_data['theta_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
        ax1.set_title('Mean potential temp')
        ax1.set_ylabel('height (m)')
        ax1.set_xlabel('potential T (K)')
        
        lgd = ax1.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
        p.savefig(testing+'_'+test_run+'theta_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        p.close(theta_fig)
        

    if testcase == 'stratus' or testcase == 'shallow_convection' :
        q_theta_fig = p.figure(figsize=(15,5))
        #q_theta_fig.suptitle(fig_title, fontsize=16 )
        ax1 = q_theta_fig.add_subplot(131)
        ax2 = q_theta_fig.add_subplot(132)
        ax3 = q_theta_fig.add_subplot(133)
        for file_index, file_numbers in enumerate(monc_dict['files']) :
            ax1.plot(monc_36Pe_data['theta_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax1.plot(monc_MultiPe_data['theta_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax2.plot(monc_36Pe_data['vapour_mmr_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax2.plot(monc_MultiPe_data['vapour_mmr_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax3.plot(monc_36Pe_data['liquid_mmr_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax3.plot(monc_MultiPe_data['liquid_mmr_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
                    
        ax1.set_title('Mean potential temp')
        ax1.set_ylabel('height (m)')
        ax1.set_xlabel('potential T (K)')
        ax2.set_title('Vapour mass mixing ratio')
        ax2.set_xlabel('qv (kg kg$^{-1}$)')
        ax3.set_title('Liquid mass mixing ratio')
        ax3.set_xlabel('ql (kg kg$^{-1}$)')

        lgd = ax3.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
        p.savefig(testing+'_'+test_run+'q_theta_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        p.close(q_theta_fig)

    if testcase == "rce" :
        q_theta_fig = p.figure(figsize=(17,10))
        #q_theta_fig.suptitle(fig_title, fontsize=16 )
        ax1 = q_theta_fig.add_subplot(231)
        ax2 = q_theta_fig.add_subplot(232)
        ax3 = q_theta_fig.add_subplot(233)
        ax4 = q_theta_fig.add_subplot(234)
        ax5 = q_theta_fig.add_subplot(235)
        ax6 = q_theta_fig.add_subplot(236)
        for file_index, file_numbers in enumerate(monc_dict['files']) :
            #ax1.plot(monc_36Pe_data['theta_mean'][0][file_index,:], z_plot[:],
            #         color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            #ax1.plot(monc_MultiPe_data['theta_mean'][0][file_index,:], z_plot[:], '--',
            #         color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax1.plot(monc_36Pe_data['vapour_mmr_mean'][0][file_index,:]*1.e3, z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax1.plot(monc_MultiPe_data['vapour_mmr_mean'][0][file_index,:]*1.e3, z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax2.plot(monc_36Pe_data['liquid_mmr_mean'][0][file_index,:]*1.e3, z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax2.plot(monc_MultiPe_data['liquid_mmr_mean'][0][file_index,:]*1.e3, z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax3.plot(monc_36Pe_data['rain_mmr_mean'][0][file_index,:]*1.e3, z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax3.plot(monc_MultiPe_data['rain_mmr_mean'][0][file_index,:]*1.e3, z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax4.plot(monc_36Pe_data['ice_mmr_mean'][0][file_index,:]*1.e3, z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax4.plot(monc_MultiPe_data[ 'ice_mmr_mean'][0][file_index,:]*1.e3, z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax5.plot(monc_36Pe_data['snow_mmr_mean'][0][file_index,:]*1.e3, z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax5.plot(monc_MultiPe_data['snow_mmr_mean'][0][file_index,:]*1.e3, z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
            ax6.plot(monc_36Pe_data['graupel_mmr_mean'][0][file_index,:]*1.e3, z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_36Pe+monc_dict['t_labels'][file_index])
            ax6.plot(monc_MultiPe_data['graupel_mmr_mean'][0][file_index,:]*1.e3, z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_72Pe+monc_dict['t_labels'][file_index])
                            
        ax1.set_title('Vapour mass mixing ratio')
        ax1.set_ylabel('height (m)')
        ax1.set_xlabel('qv (g kg$^{-1}$)')
        ax2.set_title('Liquid mass mixing ratio')
        ax2.set_xlabel('ql (g kg$^{-1}$)')
        ax3.set_title('Rain mass mixing ratio')
        ax3.set_xlabel('qr (g kg$^{-1}$)')
        ax4.set_title('Ice mass mixing ratio')
        ax4.set_ylabel('height (m)')
        ax4.set_xlabel('qi (g kg$^{-1}$)')
        ax5.set_title('Snow mass mixing ratio')
        ax5.set_xlabel('qs (g kg$^{-1}$)')
        ax6.set_title('Graupel mass mixing ratio')
        ax6.set_xlabel('qg (g kg$^{-1}$)')

        lgd = ax3.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
        p.savefig(testing+'_'+test_run+'q_mass_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        p.close(q_theta_fig) 
                     
     

    
