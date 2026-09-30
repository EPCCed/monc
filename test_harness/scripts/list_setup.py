#!/usr/bin/env python

import sys

expected_decompostion_failures = ['ScFull_2M_Ndfix_', 'ScFull_Socrates_2M_Ndfix_',
                                  'CuFull_2M_Ndfix_', 'RCE_2M_Ndfix_', 'RCESocrates_2M_Ndfix_',
                                  'ScFull_', 'ScNoRadGal_', 'CuFull_']

def main_test_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        test_list = ['WarmPw_', 'WarmPwNoSmag_', 
                     'WarmTvd_', 'WarmTvdNoSmag_',
                     'ColdPw_', 'ColdPwNoSmag_',   
                     'ColdTvd_', 'ColdTvdNoSmag_',
                     'WarmPwNoSmagGal_', 'WarmTvdNoSmagGal_', 
                     'ColdPwNoSmagGal_', 'ColdTvdNoSmagGal_',
                     'ColdNoSmagGalAdv_', 'ColdNoSmagGalMomAdv_',
                     'WarmNoSmagGalAdv_', 'WarmNoSmagGalMomAdv_' ] 

    elif case == 'stratus':
        test_list = ['ScFull_' , 'ScNoGalSrfRadSubDamp_', 'ScNoSrfSubDamp_',
                     'ScNoDamp_', 'ScNoSrfRadSubDamp_', 'ScNoSubDamp_', 'ScNoRadGal_'] 
        
    elif case == 'drybl':
        test_list = ['DryBlFull_', 
                     'DryBlNoSmagGal_', 
                     'DryBlNoGal_']

    elif case == 'shallow_convection' :
        test_list = ['CuFull_', 'CuNoDamp_', 'CuNoGalSrfForceSubDamp_', 
                     'CuNoSrfForceSubDamp_', 'CuNoSrfSubDamp_', 'CuNoSubDamp_']
        #test_list = ['CuFull_', 'CuNoDamp_']
    else :
        sys.exit(f'testcase {case} is not defined in test_names. Please check testcase name, EXITING')
    return test_list

def ecse_test_names(case):
    if case == 'bubble':
        test_list = ['WarmPw_']

    elif case == 'stratus':
        test_list = ['ScNoDamp_', 'ScNoDamp_2M_Ndfix_']

    elif case == 'shallow_convection':
        test_list = ['CuNoDamp_', 'CuNoDamp_2M_Ndfix_']

    elif case == 'rce':
        test_list = ['RCENoDampSocrates_2M_Ndfix_']
    else :
        sys.exit(f'testcase {case} is not defined in test_names. Please check testcase name, EXITING')
    return test_list

def casim_test_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus':
        test_list = ['ScFull_2M_Ndfix_', 'ScNoSubDamp_2M_Ndfix_', 'ScNoDamp_2M_Ndfix_',  
                     'ScFull_Socrates_2M_Ndfix_', 'ScNoSubDamp_Socrates_2M_Ndfix_', 'ScNoDamp_Socrates_2M_Ndfix_'] 
    elif case == 'stratus_diurnal':
        test_list = ['ScNoDamp_2M_Ndfix_diurnal_']     
    elif case == 'shallow_convection' :
        test_list = ['CuFull_2M_Ndfix_', 'CuNoDamp_2M_Ndfix_', 'CuNoDamp_Socrates_2M_Ndfix_', 'CuNoSubDamp_2M_Ndfix_']
    elif case == 'rce' :
        test_list = ['RCE_2M_Ndfix_', 'RCENoDamp_2M_Ndfix_', 'RCENoDampForce_2M_Ndfix_',
                     'RCENoDampSocrates_2M_Ndfix_', 'RCENoDampNoUVforce_2M_Ndfix_', 
                     'RCESocrates_2M_Ndfix_']
    else :
        sys.exit(f'testcase {case} is not defined in test_names. Please check testcase name, EXITING')
    return test_list

def casim_aeroproc_test_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus_hilletal':
        test_list = ['ScFull_2M_fullproc_iopt3_','ScFull_2M_noproc_iopt3_', 'ScFull_2M_passiveproc_iopt3_'] 
    elif case == 'lem_bomex' :
        test_list = ['CuNoDamp_2M_fullproc_iopt3_', 'CuNoDamp_2M_noproc_iopt3_', 'CuNoDamp_2M_passiveproc_iopt3_']     
    elif case == 'rce_deep' :
        test_list = ['RCENoDamp_2M_fullproc_iopt3_','RCENoDamp_2M_noproc_iopt3_']
    else :
        sys.exit(f'testcase {case} is not defined in test_names. Please check testcase name, EXITING')
    return test_list



def main_figure_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        fig_title_list = [ 'warm bubble - P-W mom advection, Galilean Trans & Smag mixing on',
                           'warm bubble - P-W mom advection, Galilean Trans', 
                           'warm bubble - TVD mom advection, Galilean Trans & Smag mixing on',
                           'warm bubble - TVD mom advection, Galilean Trans', 
                           'cold bubble - P-W mom advection, Galilean Trans & Smag mixing on',
                           'cold bubble - P-W mom advection, Galilean Trans',
                           'cold bubble - TVD mom advection, Galilean Trans & Smag mixing on',
                           'cold bubble - TVD mom advection, Galilean Trans',
                           'warm bubble - pressure, P-W momentum, TVD scalars',
                           'warm bubble - pressure, TVD momentum, TVD scalars',
                           'cold bubble - pressure, P-W momentum, TVD scalars',
                           'cold bubble - pressure, TVD momentum, TVD scalars',
                           'cold bubble - pressure, NO advection',
                           'cold bubble - pressure, TVD scalar advection, NO momentum advection',
                           'warm bubble - pressure, NO advection',
                           'warm bubble - pressure, TVD scalar advection, NO momentum advection',
                           'warm bubble - pressure, NO scalar advection, PW momentum advection',
                           'cold bubble - pressure, NO scalar advection, PW momentum advection',
                           'cold bubble - pressure, PW scalar and momentum advection',
                           'warm bubble - pressure, PW scalar and momentum advection']
    elif case == 'stratus':
        fig_title_list = [ 'Stratus - Full',
                           'Stratus - No Galilean, surface flux, radiation param, subsidence, damping',
                           'Stratus - No surface flux, subsidence, damping', 
                           'Stratus - No damping',
                           'Stratus - No surface flux, radiation param, subsidence, damping',
                           'Stratus - No subsidence, damping', 
                           'Stratus - No radiation, galilean transformation']
            
    elif case == 'drybl':
        fig_title_list = [ 'Dry Boundary layer - Full',
                           'Dry Boundary layer - No Smag or Galilean trans',
                           'Dry Boundary layer - No Galilean']
    elif case == 'shallow_convection':
        fig_title_list = [ 'Shallow convection - Full',
                           'Shallow convection - No damping',
                           'Shallow convection - No Galilean, surface flux, forcing, subsidence, damping', 
                           'Shallow convection - No surface flux, forcing, subsidence, damping',
                           'Shallow convection - No surface flux, subsidence, damping',
                           'Shallow convection - No subsidence, damping']
    else :
        sys.exit(f'testcase {case} is not defined in figure_names. Please check testcase name, EXITING')
    return fig_title_list

def ecse_figure_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        fig_title_list = [ 'warm bubble - P-W mom advection, Galilean Trans & Smag mixing on']
    elif case == 'stratus':
        fig_title_list = [ 'Stratus - No damping',
                           'Stratus - No damping, casim, Nd = 50']
    elif case == 'shallow_convection':
        fig_title_list = [ 'Shallow convection - No damping',
                           'Shallow convection - No damping, casim, Nd = 50']
    elif case == 'rce':
        fig_title_list = [ 'RCE - No Damp, socrates, casim, Nd = 50']
    else :
        sys.exit(f'testcase {case} is not defined in figure_names. Please check testcase name, EXITING')
    return fig_title_list


def casim_figure_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus':
        fig_title_list = [ 'Stratus - Full, casim, Nd = 50',
                           'Stratus - No subsidence, damping, casim, Nd = 50', 
                           'Stratus - No damping, casim, Nd = 50', 
                           'Stratus - Full, casim, socreate, Nd = 50',
                           'Stratus - No subsidence, damping, casim, socrates, Nd = 50', 
                           'Stratus - No damping, casim, socrates, Nd = 50']
    elif case == 'stratus_diurnal' :
        fig_title_list = ['Diurnal Stratus - No damping, casim, Nd = 50'] 
    elif case == 'shallow_convection':
        fig_title_list = [ 'Shallow convection - Full, casim, Nd = 50',
                           'Shallow convection - No damping, casim, Nd = 50', 
                           'Shallow convection - No damping, casim, socrates, Nd = 50',
                           'Shallow convection - No subsidence, damping, casim, Nd = 50']
    elif case == 'rce' :
        fig_title_list = [ 'RCE - Full, prescribed cooling, casim, Nd = 50',
                           'RCE - No Damp , casim, Nd = 50', 
                           'RCE - No Damp, no force, casim, Nd = 50',
                           'RCE - No Damp, socrates, casim, Nd = 50',
                           'RCE - No Damp, no wind relax, socrates, casim, Nd = 50', 
                           'RCE - Full, socrates, casim, Nd = 50']
    else :
        sys.exit(f'testcase {case} is not defined in figure_names. Please check testcase name, EXITING')
    return fig_title_list

def casim_aeroproc_figure_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus_hilletal':
        fig_title_list = [ 'Stratus - Na = 100, ARG, full proc',
                           'Stratus - Na = 100, ARG, no proc', 
                           'Stratus - Na = 100, ARG, passive proc'] 
    elif case == 'lem_bomex':
        fig_title_list = [ 'Shallow convection - Na = 100, ARG, full proc'
                           'Shallow convection - Na = 100, ARG, no proc', 
                           'Shallow convection - Na = 100, ARG, passive proc']
    elif case == 'rce' :
        fig_title_list = [ 'RCE - Na = 100, ARG, full proc',
                           'RCE - Na = 100, ARG, no proc'] 
    else :
        sys.exit(f'testcase {case} is not defined in figure_names. Please check testcase name, EXITING')
    return fig_title_list
