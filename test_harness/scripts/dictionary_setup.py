#!/usr/bin/env python

# Used by: monc_kgo_bit_compare.py
#          monc_speed_test.py
#          monc_vs_lem.py

import sys

def monc_main(case, kgo_datadir, new_datadir) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_3600.nc', 'diagnostics_3600.nc', 'diagnostics_3600.nc', 
                          'diagnostics_3600.nc', 'diagnostics_3600.nc', 'diagnostics_3600.nc'],
                     't_labels':['10 mins', '20 mins', '30 mins','40 mins', '50 mins', '60 mins'],
                     'time_color':[ 'c', 'm', 'k', 'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    elif case == 'stratus':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_7260.nc', 'diagnostics_7260.nc', 'diagnostics_7260.nc',
                              'diagnostics_7260.nc'],
                     't_labels':['30 mins', '60 mins', '90 mins', '120 min'],
                     'time_color':[ 'blue', 'red', 'y', 'c' ],
                     'no_procs':['36_', '72_']}            
    elif case == 'drybl':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_32400.nc', 
                              'diagnostics_32400.nc', 'diagnostics_32400.nc', 'diagnostics_32400.nc', 
                              'diagnostics_32400.nc', 'diagnostics_32400.nc', 'diagnostics_32400.nc'],
                     't_labels':[' 3 hours', ' 4 hours',
                                 ' 5 hour', ' 6 hours', ' 7 hours', ' 8 hours',
                                 ' 9 hours'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y', 'c', ],
                     'no_procs':['36_', '72_']}
    elif case == 'shallow_convection':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc',
                              'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins', 
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    else :
        sys.exit('testcase is not defined. Please check testcase name, EXITING')
    return monc_dict

def monc_casim(case, kgo_datadir, new_datadir) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_7260.nc', 'diagnostics_7260.nc', 'diagnostics_7260.nc',
                              'diagnostics_7260.nc'],
                     't_labels':['30 mins', '60 mins', '90 mins', '120 min'],
                     'time_color':[ 'blue', 'red', 'y', 'c' ],
                     'no_procs':['36_', '72_']} 
    elif case == 'stratus_diurnal':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_86400.nc', 'diagnostics_86400.nc', 'diagnostics_86400.nc',
                              'diagnostics_86400.nc'],
                     't_labels':['22.5 hours', '23 hours', '23.5 hours', '24 hours'],
                     'time_color':[ 'blue', 'red', 'y', 'c' ],
                     'no_procs':['36_', '72_']} 
    elif case == 'shallow_convection':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc',
                              'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins', 
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    elif case == 'rce':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_86400.nc', 'diagnostics_86400.nc', 'diagnostics_86400.nc'],
                     't_labels':['22 hours', '23 hours', '24 hours'],
                     'time_color':[ 'm', 'g', 'k'],
                     'no_procs':['36_', '72_']}   
    else :
        sys.exit('testcase is not defined. Please check testcase name, EXITING')
    return monc_dict

def monc_casim_aeroproc(case, kgo_datadir, new_datadir) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus_hilletal':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_7260.nc', 'diagnostics_7260.nc', 'diagnostics_7260.nc',
                              'diagnostics_7260.nc'],
                     't_labels':['30 mins', '60 mins', '90 mins', '120 min'],
                     'time_color':[ 'blue', 'red', 'y', 'c' ],
                     'no_procs':['36_', '72_']} 
    elif case == 'lem_bomex':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc',
                              'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins', 
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    elif case == 'rce_deep':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_86400.nc', 'diagnostics_86400.nc', 'diagnostics_86400.nc'],
                     't_labels':['22 hours', '23 hours', '24 hours'],
                     'time_color':[ 'm', 'g', 'k'],
                     'no_procs':['36_', '72_']}   
    else :
        sys.exit('testcase is not defined. Please check testcase name, EXITING')
    return monc_dict

def lem(case, kgo_datadir):
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        lem_dict ={'kgo_dir' :kgo_datadir,
                   'files':['diagnostics_0002', 'diagnostics_0003', 'diagnostics_0004', 
                            'diagnostics_0005', 'diagnostics_0005', 'diagnostics_0007' ], 
                   't_labels':['LEM 10 mins', 'LEM 20 mins', 'LEM 30 mins','LEM 40 mins', 
                               'LEM 50 mins', 'LEM 60 mins' ],
                   'no_procs':['2_', '36_']}
    elif case == 'stratus':
        lem_dict ={'kgo_dir' :kgo_datadir,
                   'files':['diagnostics_0002', 'diagnostics_0003', 'diagnostics_0004', 
                            'diagnostics_0005'], 
                   't_labels':['LEM 30 mins', 'LEM 60 mins', 'LEM 90 mins', 'LEM 120 min'],
                   'no_procs':['2_', '32_']}
        
    elif case == 'drybl':
        lem_dict ={'kgo_dir' :kgo_datadir,
                   'files':['diagnostics_0004', 
                            'diagnostics_0005', 'diagnostics_0006', 'diagnostics_0007', 
                            'diagnostics_0008', 'diagnostics_0009', 'diagnostics_0010'], 
                   't_labels':['LEM 3 hours', 'LEM 4 hours',
                               'LEM 5 hour', 'LEM 6 hours', 'LEM 7 hours', 'LEM 8 hours',
                               'LEM 9 hours'],
                   'no_procs':['2_', '32_']}
    
    else :
        sys.exit('testcase is not defined in dictionary_setup.py. Please check testcase name, EXITING')
    return lem_dict

def monc_ecse(case, kgo_datadir, new_datadir):
    # call module that sets up the lists and dictionaries for the case    
    if case == 'bubble':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_3600.nc', 'diagnostics_3600.nc', 'diagnostics_3600.nc',
                          'diagnostics_3600.nc', 'diagnostics_3600.nc', 'diagnostics_3600.nc'],
                     't_labels':['10 mins', '20 mins', '30 mins','40 mins', '50 mins', '60 mins'],
                     'time_color':[ 'c', 'm', 'k', 'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    elif case == 'stratus':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_7260.nc', 'diagnostics_7260.nc', 'diagnostics_7260.nc',
                              'diagnostics_7260.nc'],
                     't_labels':['30 mins', '60 mins', '90 mins', '120 min'],
                     'time_color':[ 'blue', 'red', 'y', 'c' ],
                     'no_procs':['36_', '72_']}
    elif case == 'shallow_convection':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc',
                              'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins',
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    elif case == 'rce':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins',
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['36_', '72_']}
    else :
        sys.exit('testcase is not defined. Please check testcase name, EXITING')
    return monc_dict

