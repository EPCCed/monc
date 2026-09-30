#!/usr/bin/env python

# Beta test-suite (see ticket #65)
# Code for extracting data and plotting the comparison of the LEM and MONC 
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
# Adhill - 040516

import matplotlib.pyplot as p
import sys
import os
#import local module
import get_data

compiler = 'cray'
kgo_dir = 'latest_kgo/r5815/'
suite = 'main_component'

# read in lists and dictionaries, which are
# required to setup and run the test_harness
dict_list = test_harness_setup.dict_list_read(suite)
procs = dict_list[0]
test_case = dict_list[1]

new_monc_datadir = input("Please provide full path for directory with test data: ")

if not os.path.isdir(new_monc_datadir) :
    sys.exit(new_monc_datadir+' does not exist, please check and try again - Exiting')
    
######### Qualitatively compare LEM and MONC output from component testing #############################

if optimisation_level == 1 :
    monc_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt1/monc/'
else :
    monc_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt3/monc/'+kgo_dir+'/'
    
lem_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/cray/opt1/lem/'

for case in test_case :
    
    test_list = list_setup.test_names(case)
    figure_title_list = list_setup.figure_names(case)
    
    monc_dict = dictionary_setup.monc(case, monc_kgo_datadir, new_monc_datadir)
    lem_dict = dictionary_setup.lem(case, lem_kgo_datadir)
        
    monc_vs_lem.plotting(case, test_list, figure_title_list, monc_dict, lem_dict)

def plotting(testcase, test_list, fig_title_list, monc_dict, lem_dict):
    
    monc_proc_numbers = '36_'
    lem_proc_numbers = '32_'
    linestyle_for_time='o'

    for test_idx, test_run in enumerate(test_list) :
        # first extract all the data required and put into a dictionary

        model_fn = monc_dict['new_dir']+test_run+monc_proc_numbers
        # this dictionary is used to output data from get_profile 
        #it is overwritten on each loop
        monc_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                      'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                      'vapour_mmr_mean': [], 'liquid_mmr_mean': []}
        print(model_fn)
        get_data.profile_diagnostics(model_fn, monc_dict['files'], monc_data)

        model_fn = lem_dict['kgo_dir']+test_run+lem_proc_numbers                
        # this dictionary is used to output data from get_profile. 
        # It overwritten on each loop
        lem_data = { 'z':[], 'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                     'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],  
                     'vapour_mmr_mean': [], 'liquid_mmr_mean': [] }

        get_data.profile_diagnostics(model_fn, lem_dict['files'], lem_data) 
        # Assumes the z dimension uses the z from the LEM for both LEM and MONC    
        z_plot = lem_data['z'][0][:]

        # Set up the velocity figure for the case 
        #
        velocity_fig = p.figure(figsize=(15,10))
        velocity_fig.suptitle(fig_title_list[test_idx], fontsize=16 )
        ax1 = velocity_fig.add_subplot(231)
        ax2 = velocity_fig.add_subplot(232)
        ax3 = velocity_fig.add_subplot(233)
        ax4 = velocity_fig.add_subplot(234)
        ax5 = velocity_fig.add_subplot(235)
        ax6 = velocity_fig.add_subplot(236)

        for file_index, file_numbers in enumerate(monc_dict['files']) :
            if testcase == 'drybl' :
                monc_index = file_index+2
                lem_index = file_index
            else :
                monc_index = file_index
                lem_index = file_index
            
            print(file_index)
            ax1.plot(monc_data['u_mean'][0][monc_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index],label='MONC '+monc_dict['t_labels'][file_index])
            ax1.plot(lem_data['u_mean'][file_index][:], z_plot[:]
                     , '--',color=monc_dict['time_color'][file_index],label='LEM '+monc_dict['t_labels'][file_index] )
            
            ax2.plot(monc_data['v_mean'][0][monc_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
            ax2.plot(lem_data['v_mean'][file_index][:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
            
            ax3.plot(monc_data['w_mean'][0][monc_index,:], z_plot[:],
                    color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
            ax3.plot(lem_data['w_mean'][file_index][:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
            
            ax4.plot(monc_data['uu_mean'][0][monc_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
            ax4.plot(lem_data['uu_mean'][file_index][:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
            
            ax5.plot(monc_data['vv_mean'][0][monc_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
            ax5.plot(lem_data['vv_mean'][file_index][:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
                        
            ax6.plot(monc_data['ww_mean'][0][monc_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
            ax6.plot(lem_data['ww_mean'][file_index][:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
    
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
        
        p.savefig(test_run+'velocity_monc_vs_lem.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

        p.close(velocity_fig)

        if testcase == 'bubble' :
            theta_fig = p.figure(figsize=(8,5))
            theta_fig.suptitle(fig_title_list[test_idx], fontsize=16 )
            ax1 = theta_fig.add_subplot(111)
            for file_index, file_numbers in enumerate(monc_dict['files']) :
                ax1.plot(monc_data['theta_mean'][0][file_index,:], z_plot[:],
                         color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
                ax1.plot(lem_data['theta_mean'][file_index][:], z_plot[:], '--',
                         color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
            ax1.set_title('Mean potential temp')
            ax1.set_ylabel('height (m)')
            ax1.set_xlabel('potential T (K)')
            lgd = ax1.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
            p.savefig(test_run+'theta_monc_vs_lem.png', bbox_extra_artists=(lgd,), bbox_inches='tight')

            p.close(theta_fig)

        if testcase == 'stratus' :
            q_theta_fig = p.figure(figsize=(15,5))
            q_theta_fig.suptitle(fig_title_list[test_idx], fontsize=16 )
            ax1 = q_theta_fig.add_subplot(131)
            ax2 = q_theta_fig.add_subplot(132)
            ax3 = q_theta_fig.add_subplot(133)
            for file_index, file_numbers in enumerate(monc_dict['files']) :
                ax1.plot(monc_data['theta_mean'][0][file_index,:], z_plot[:],
                         color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
                ax1.plot(lem_data['theta_mean'][file_index][:], z_plot[:], '--',
                         color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
                ax2.plot(monc_data['vapour_mmr_mean'][0][file_index,:], z_plot[:],
                         color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
                ax2.plot(lem_data['vapour_mmr_mean'][file_index][:], z_plot[:], '--',
                         color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
                ax3.plot(monc_data['liquid_mmr_mean'][0][file_index,:], z_plot[:],
                         color=monc_dict['time_color'][file_index], label='MONC '+monc_dict['t_labels'][file_index])
                ax3.plot(lem_data['liquid_mmr_mean'][file_index][:], z_plot[:], '--',
                         color=monc_dict['time_color'][file_index], label='LEM '+monc_dict['t_labels'][file_index])
                    
            ax1.set_title('Mean potential temp')
            ax1.set_ylabel('height (m)')
            ax1.set_xlabel('potential T (K)')
            ax2.set_title('Vapour mass mixing ratio')
            ax2.set_xlabel('qv (kg kg$^{-1}$')
            ax3.set_title('Liquid mass mixing ratio')
            ax3.set_xlabel('ql (kg kg$^{-1}$')

            lgd = ax3.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
            p.savefig(test_run+'q_theta_monc_vs_lem.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        
            p.close(q_theta_fig)

# 
            
