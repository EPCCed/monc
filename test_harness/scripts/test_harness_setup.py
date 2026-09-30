#!/usr/bin/env python

import sys

def dict_list_read(suite) :
    # call module that sets up the lists and dictionaries for the case
    print('\nSetting up:')
    if suite == 'main_component' :
        print( 'main_component')
        procs = [72, 36] 
        testcases = [ 'bubble' , 'drybl', 'stratus', 'shallow_convection' ]
        
        test_cfg={ 'bubble':[ 'ColdNoSmagGalAdv', 'ColdNoSmagGalMomAdv',
                              'ColdPw','ColdPwNoSmag','ColdPwNoSmagGal', 
                              'ColdTvd','ColdTvdNoSmag','ColdTvdNoSmagGal',
                              'WarmNoSmagGalAdv', 'WarmNoSmagGalMomAdv',
                              'WarmPw', 'WarmPwNoSmag', 'WarmPwNoSmagGal',  
                              'WarmTvd', 'WarmTvdNoSmag', 'WarmTvdNoSmagGal'],
                   'drybl':['DryBlFull', 'DryBlNoSmagGal', 'DryBlNoGal'],
                   'stratus':['ScFull', 'ScNoDamp', 'ScNoSrfRadSubDamp', 'ScNoRadGal', 
                              'ScNoGalSrfRadSubDamp', 'ScNoSrfSubDamp', 'ScNoSubDamp', 
                              'ScNoRad', 'ScFixFluxNoSubDamp'], 
                   'shallow_convection':[ 'CuFull', 'CuNoDamp' , 'CuNoGalSrfForceSubDamp', 
                                          'CuNoSrfForceSubDamp', 'CuNoSrfSubDamp', 'CuNoSubDamp'] 
        }
    elif suite == 'standard' :
        procs = [36]
        testcases = ['tank_experiments', 'drybl', 'stratus', 'shallow_convection', 'stable',
                     'radiative_convective_equilibrium']
        test_cfg={'tank_experiments':['cold_bubble', 'warm_bubble'],
                  'drybl':['drybl'],
                  'stratus':['mbl_sc', 'mbl_sc_casim', 'mbl_sc_diurnal'],
                  'shallow_convection':['bomex', 'bomex_casim', 'bomex_casim_socrates'],
                  'stable':['lanfex_IOP1_casim'],
                  'radiative_convective_equilibrium':['RCE_casim', 'RCE_casim_socrates']
        }
    elif suite == 'casim_socrates':
        print( 'casim_socrates')
        procs = [72,36] 
        testcases = [ 'stratus', 'stratus_diurnal', 'shallow_convection', 'rce' ]
        
        test_cfg={ 
            'stratus':[ 'ScNoSubDamp_2M_Ndfix', 'ScNoDamp_2M_Ndfix', 'ScFull_2M_Ndfix', 
                        'ScNoSubDamp_Socrates_2M_Ndfix', 'ScNoDamp_Socrates_2M_Ndfix', 'ScFull_Socrates_2M_Ndfix'],
            'stratus_diurnal':['ScNoDamp_2M_Ndfix_diurnal'],
            'shallow_convection':['CuNoSubDamp_2M_Ndfix','CuNoDamp_2M_Ndfix', 'CuNoDamp_Socrates_2M_Ndfix', 'CuFull_2M_Ndfix'],
            'rce':['RCENoDampSocrates_2M_Ndfix','RCENoDampForce_2M_Ndfix', 'RCENoDamp_2M_Ndfix', 'RCE_2M_Ndfix', 
                   'RCESocrates_2M_Ndfix', 'RCENoDampNoUVforce_2M_Ndfix']
        }
    elif suite == 'casim_aerosol_processing' :
        print( 'casim_socrates_processing')
        procs = [72, 36] 
        testcases = [ 'stratus_hilletal', 'lem_bomex', 'rce_deep' ]
        
        test_cfg={ 
            'stratus_hilletal':[ 'ScFull_2M_fullproc_iopt3',  'ScFull_2M_fullproc_iopt5',  'ScFull_2M_NdFix',
                                 'ScFull_2M_noproc_iopt3',  'ScFull_2M_noproc_iopt5',  'ScFull_2M_passiveproc_iopt3',
                                 'ScFull_2M_passiveproc_iopt5' ],
            'lem_bomex':['CuNoDamp_2M_fullproc_iopt3',  'CuNoDamp_2M_fullproc_iopt5',  'CuNoDamp_2M_NdFix',
                         'CuNoDamp_2M_noproc_iopt3',  'CuNoDamp_2M_noproc_iopt5',  'CuNoDamp_2M_passiveproc_iopt3',
                         'CuNoDamp_2M_passiveproc_iopt5'],
            'rce_deep':['RCENoDamp_2M_fullproc_iopt3',  'RCENoDamp_2M_fullproc_iopt5',  
                        'RCENoDamp_2M_noproc_iopt3',  'RCENoDamp_2M_noproc_iopt5',  'RCENoDamp_2M_passiveproc_iopt3',
                        'RCENoDamp_2M_passiveproc_iopt5']
        }
    elif suite == 'ecse':
        print('ecse')
        procs = [72, 36]
        testcases = [ 'bubble', 'stratus', 'shallow_convection', 'rce' ]

        test_cfg={
            'bubble':
                    ['WarmPw'],
            'stratus':
                    ['ScNoDamp', 'ScNoDamp_2M_Ndfix'],
            'shallow_convection':
                    ['CuNoDamp', 'CuNoDamp_2M_Ndfix'],
            'rce':
                    ['RCENoDampSocrates_2M_Ndfix']
        }
    else :
        sys.exit(f"\n Suite '{suite}' is not defined. Please check suite name, EXITING")

    return  procs, testcases, test_cfg        
