import os,shutil

def gen_fcdft_input(xyz,fc_param,parameterization,charges,xc,bs,scf_threshold):
# Generate the FCDFT input files from a list of arguments
# 'xyz', type=str, default=[None], nargs='*', help='molecule xyz files, eg 'h2.xyz''
# 'fc_param', type=list, default=[None], help='list of fc parameters'
# 'parameterization', type=str, default='calc', help='parameterization or use given parameters'
# 'charge', type=int, default=[], nargs='*', help='charge'
# 'xc', type=str, default='PBE', help='DFT functional'
# 'bs', type=str, default='STO-3G', help='basis set'

    molecules = []
    # name of the molecule(s), eg 'h2'
    for f in xyz:
        molecules.append(f.split('.')[0])

    if charges == []:
        charges = [0]*len(molecules)

    for i in range(0,len(molecules)):
        genInput(molecules[i], charges[i], bs, xc, parameterization, scf_threshold, fc_param)


def genInput(molecule,charge,bs,xc,parameterization,scf_threshold,fc_param):
    current_dir = os.getcwd()
    dirName = current_dir + '/' + 'FC' + xc + '_' + parameterization
    xyz_dir = current_dir + '/xyz'

    if not os.path.exists(dirName):
        os.makedirs(dirName)
    shutil.copy2(xyz_dir + '/' + molecule+'.xyz', dirName + '/' + molecule + '.xyz')

    # use Kaito's pairwise potentials
    if parameterization == 'PBE':
        if xc == 'LDA':
            # Kaito: training={H,C}, LDA/STO-3G_PBE/6-31G*
            shutil.copy2(current_dir + '/' + 'u_fclda_pbe.txt', dirName + '/' + 'u_fclda_pbe.txt')
        if xc == 'Hartree':
            # Kaito: training={H,C}, Hartree/STO-3G_PBE/6-31G*
            shutil.copy2(current_dir + '/' + 'u_fch_pbe.txt', dirName + '/' + 'u_fch_pbe.txt')
    elif parameterization == 'B3LYP':
        if xc == 'LDA':
            # Kaito: training={H,C}, LDA/STO-3G_B3LYP/6-31G*
            shutil.copy2(current_dir + '/' + 'u_fclda_b3lyp.txt', dirName + '/' + 'u_fclda_b3lyp.txt')
        if xc == 'Hartree':
            # Kaito: training={H,C}, Hartree/STO-3G_B3LYP/6-31G*
            shutil.copy2(current_dir + '/' + 'u_fch_b3lyp.txt', dirName + '/' + 'u_fch_b3lyp.txt')

    with open(os.path.join(dirName,molecule+'.in'),'w') as f:
        f.write(dft % (molecule,charge,xc,bs,scf_threshold))
        # parameterization
        if parameterization == 'calc':
            # shutil.copy2(current_dir + '/' + 'u_fclda_pbe.txt', dirName + '/' + 'u_fclda_pbe.txt')

            # H only
            # f.write(fc_variable_h % (str(fc_param[0]),str(fc_param[1]),str(fc_param[2]))) 
            # {H,C}
            f.write(fc_variable % (\
                str(fc_param[0]),str(fc_param[1]),str(fc_param[2]),str(fc_param[3]),\
                str(fc_param[4]),str(fc_param[5]),str(fc_param[6]) ) )

        elif parameterization == 'Ucalc':
            # H only with U optimized simutaneously
            f.write(fc_variable_u_h % (str(fc_param[0]),str(fc_param[1]),str(fc_param[2]))) 
            with open("FCLDA_Ucalc/u_H.txt", "w") as u_file:
                #write u file
                u_file.write(U_variable_h % (\
                        str(fc_param[3]),str(fc_param[4]),str(fc_param[5]),str(fc_param[6]),str(fc_param[7]),\
                        str(fc_param[8]),str(fc_param[9]),str(fc_param[10]),str(fc_param[11]),str(fc_param[12]),\
                        str(fc_param[13]),str(fc_param[14]),str(fc_param[15]),str(fc_param[16]),str(fc_param[17]),\
                        str(fc_param[18]),str(fc_param[19]),str(fc_param[20]),str(fc_param[21]),str(fc_param[22]),))

        elif parameterization == 'U_only_calc':
            # H only with U optimized simutaneously
            f.write(fc_variable_u_h % ('-0.0408557613', '0.0520289523', '0.9850066033')) 
            with open("FCLDA_U_only_calc/u_H.txt", "w") as u_file:
                #write u file
                u_file.write(U_variable_h % (\
                        str(fc_param[0]),str(fc_param[1]),str(fc_param[2]),\
                        str(fc_param[3]),str(fc_param[4]),str(fc_param[5]),str(fc_param[6]),str(fc_param[7]),\
                        str(fc_param[8]),str(fc_param[9]),str(fc_param[10]),str(fc_param[11]),str(fc_param[12]),\
                        str(fc_param[13]),str(fc_param[14]),str(fc_param[15]),str(fc_param[16]),str(fc_param[17]),\
                        str(fc_param[18]),str(fc_param[19])))
        
        elif parameterization == 'calc_diag':
            # H only with U optimized simutaneously
            f.write(fc_variable_diag % (str(fc_param[0]),str(fc_param[1]),str(fc_param[2]),str(fc_param[3])))

        # use existing parameters
        elif parameterization == 'PBE':
            if xc == 'LDA':
                f.write(fc_pbe_lda)
            elif xc == 'Hartree':
                f.write(fc_pbe_hartree)
            else:
                print('Unrecognized exchange-correlation functional!')
        elif parameterization == 'B3LYP':
            if xc == 'LDA':
                f.write(fc_b3lyp_lda)
            elif xc == 'Hartree':
                f.write(fc_b3lyp_hartree)
            else:
                print('Unrecognized exchange-correlation functional!')
        elif parameterization == 'diag':
            f.write(fc_sto3g_pbe_631gd_old)
        else:
            print('Unrecognized parameterization option!')
        f.write(tail)

 
# dft options in entos
#        density_fitting = '%s' true or false. 
#        print_level = %d different output options. larger->more information printed

dft ="""dft(
    structure( xyz = '%s.xyz')
    charge = '%s'
    xc = '%s'
    basis = '%s'
    print_level = 1
    energy_threshold = %s
"""
#pairwise_potential = 'u_fclda_pbe.txt'
fc_variable_h = """    
    fock_correction(
        diagonal = {'H.1s' : %s
                    }
        alpha = {'H.1s.H.1s' : %s
                 }
        beta = {'H.1s.H.1s' : %s
                }
        gamma = {'H.1s.H.1s' : 2.0
                }
        delta = {
                }
        )
"""
U_variable_h = """! Parameters for pairwise potential with 4th-order spline interpolation
! FCLDA/STO-3G fitted to PBE/6-31G*

****
H H
5 5
  1.300 %s            %s            %s            %s            %s          
  1.400 %s            %s            %s            %s            %s          
  1.600 %s            %s            %s            %s            %s          
  1.800 %s            %s            %s            %s            %s          
  2.000 0.0000000000  0.0000000000  0.0000000000  0.0000000000  0.0000000000
****
"""
fc_variable_u_h = """    pairwise_potential = 'u_H.txt' 
    fock_correction(
        diagonal = {'H.1s' : %s
                    }
        alpha = {'H.1s.H.1s' : %s
                 }
        beta = {'H.1s.H.1s' : %s
                }
        gamma = {'H.1s.H.1s' : 2.0
                }
        delta = {
                }
        )
"""
#pairwise_potential = 'u_fclda_pbe.txt'
fc_variable = """    
    fock_correction(
        diagonal = {'H.1s' : %s,
                    'C.1s' : %s,
                    'C.2s' : %s,
                    'C.2p' : %s
                    }
        off_diagonal = {'H.1s' : %s,
                    'C.2s' : %s,
                    'C.2p' : %s
                    }
        factor =  -1.0
        )
"""
#pairwise_potential = 'u_fclda_pbe.txt'
fc_variable_hc = """    pairwise_potential = 'u_fclda_pbe.txt'
    fock_correction(    
        diagonal = {'H.1s' : -0.084080685,
                    'C.1s' : %s,
                    'C.2s' : %s,
                    'C.2p' : %s
                    }
        alpha = {'H.1s.H.1s' : 0.17584356,
                 'H.1s.C.1s' : %s,
                 'H.1s.C.2s' : %s,
                 'H.1s.C.2p' : %s,
                 'C.1s.C.1s' : %s,
                 'C.1s.C.2s' : %s,
                 'C.1s.C.2p' : %s,
                 'C.2s.C.2s' : %s,
                 'C.2s.C.2p' : %s,
                 'C.2p.C.2p' : %s,
                 'C.2p.C.2p.pi' : %s
                   }
        beta = {'H.1s.H.1s' : 0.968401976,
                'H.1s.C.1s' : %s,
                'H.1s.C.2s' : %s,
                'H.1s.C.2p' : %s,
                'C.1s.C.1s' : %s,
                'C.1s.C.2s' : %s,
                'C.1s.C.2p' : %s,
                'C.2s.C.2s' : %s,
                'C.2s.C.2p' : %s,
                'C.2p.C.2p' : %s,
                'C.2p.C.2p.pi' : %s
                }
        gamma = {'H.1s.H.1s' : 2.0000000000,
                 'H.1s.C.1s' : 2.0000000000,
                 'H.1s.C.2s' : 2.0000000000,
                 'H.1s.C.2p' : 2.0000000000,
                 'C.1s.C.1s' : 2.0000000000,
                 'C.1s.C.2s' : 2.0000000000,
                 'C.1s.C.2p' : 2.0000000000,
                 'C.2s.C.2s' : 2.0000000000,
                 'C.2s.C.2p' : 2.0000000000,
                 'C.2p.C.2p' : 2.0000000000,
                 'C.2p.C.2p.pi' : 2.0000000000
                }
        delta = {'C.2p.C.2p' : %s}
        )
"""
fc_variable_diag = """    fock_correction(
        diagonal = {'H.1s' : %s,
                    'C.1s' : %s,
                    'C.2s' : %s,
                    'C.2p' : %s
                    }
        )
"""
fc_pbe_lda_brooke = """    pairwise_potential = 'u_fclda_pbe.txt'
    fock_correction(
        diagonal = {'H.1s' : %s,
                    'C.1s' : %s,
                    'C.2s' : %s,
                    'C.2p' : %s
                    }
        alpha = {'H.1s.H.1s' : %s,
                 'H.1s.C.1s' : %s,
                 'H.1s.C.2s' : %s,
                 'H.1s.C.2p' : %s,
                 'C.1s.C.1s' : %s,
                 'C.1s.C.2s' : %s,
                 'C.1s.C.2p' : %s,
                 'C.2s.C.2s' : %s,
                 'C.2s.C.2p' : %s,
                 'C.2p.C.2p' : %s,
                 'C.2p.C.2p.pi' : %s
                 }
        beta = {'H.1s.H.1s' : %s,
                'H.1s.C.1s' : %s,
                'H.1s.C.2s' : %s,
                'H.1s.C.2p' : %s,
                'C.1s.C.1s' : %s,
                'C.1s.C.2s' : %s,
                'C.1s.C.2p' : %s,
                'C.2s.C.2s' : %s,
                'C.2s.C.2p' : %s,
                'C.2p.C.2p' : %s,
                'C.2p.C.2p.pi' : %s
                }
        gamma = {'H.1s.H.1s' : 2.0000000000,
                 'H.1s.C.1s' : 2.0000000000,
                 'H.1s.C.2s' : 2.0000000000,
                 'H.1s.C.2p' : 2.0000000000,
                 'C.1s.C.1s' : 2.0000000000,
                 'C.1s.C.2s' : 2.0000000000,
                 'C.1s.C.2p' : 2.0000000000,
                 'C.2s.C.2s' : 2.0000000000,
                 'C.2s.C.2p' : 2.0000000000,
                 'C.2p.C.2p' : 2.0000000000,
                 'C.2p.C.2p.pi' : 2.0000000000
                }
        delta = {'C.2p.C.2p' : %s}
        )
"""
fc_sto3g_pbe_631gd_old = """    fock_correction(
        diagonal = {'H.1s' : -0.0579298305937,
                    'C.1s' : -0.895617860555,
                    'C.2s' : 0.0103429150385,
                    'C.2p' : -0.0743785605307,
                    'O.1s' : 0.790356507631,
                    'O.2s' : -0.234135316062,
                    'O.2p' : -0.213328780723
                    }
        )
"""
fc_sto3g_pbe_631gd = """    fock_correction(
        diagonal = {'H.1s' : -0.057695541,
                    'C.1s' : -0.739450144,
                    'C.2s' : 0.003819713,
                    'C.2p' : -0.067856055
                    }
        )
"""
fc_pbe_lda = """    pairwise_potential = 'u_fclda_pbe.txt'
    fock_correction(
        diagonal = {'H.1s' : -0.0841356082,
                    'C.1s' : -0.4514212430,
                    'C.2s' : -0.0474520673,
                    'C.2p' : -0.0968652872
                    }
        alpha = {'H.1s.H.1s' : 0.1754880254,
                 'H.1s.C.1s' : 0.0670257256,
                 'H.1s.C.2s' : -0.1651256150,
                 'H.1s.C.2p' : -0.7093263613,
                 'C.1s.C.1s' : -0.0055168832,
                 'C.1s.C.2s' : 0.0425207030,
                 'C.1s.C.2p' : -0.0182564179,
                 'C.2s.C.2s' : -0.0600625243,
                 'C.2s.C.2p' : -0.0096123275,
                 'C.2p.C.2p' : 0.0629574539,
                 'C.2p.C.2p.pi' : -0.1994883855
                 }
        beta = {'H.1s.H.1s' : 0.9695295991,
                'H.1s.C.1s' : 1.0070109726,
                'H.1s.C.2s' : 1.0000833609,
                'H.1s.C.2p' : 0.9126804570,
                'C.1s.C.1s' : 0.9999955182,
                'C.1s.C.2s' : 0.9981098167,
                'C.1s.C.2p' : 1.0134896631,
                'C.2s.C.2s' : 0.9975383311,
                'C.2s.C.2p' : 1.0016496033,
                'C.2p.C.2p' : 0.9346390221,
                'C.2p.C.2p.pi' : 0.9878326355
                }
        gamma = {'H.1s.H.1s' : 2.0000000000,
                 'H.1s.C.1s' : 2.0000000000,
                 'H.1s.C.2s' : 2.0000000000,
                 'H.1s.C.2p' : 2.0000000000,
                 'C.1s.C.1s' : 2.0000000000,
                 'C.1s.C.2s' : 2.0000000000,
                 'C.1s.C.2p' : 2.0000000000,
                 'C.2s.C.2s' : 2.0000000000,
                 'C.2s.C.2p' : 2.0000000000,
                 'C.2p.C.2p' : 2.0000000000,
                 'C.2p.C.2p.pi' : 2.0000000000
                }
        delta = {'C.2p.C.2p' : -0.252557975}
        )
"""
fc_b3lyp_lda = """    pairwise_potential = 'u_fclda_b3lyp.txt'
    fock_correction(
        diagonal = {'H.1s' : -0.1124596170,
                    'C.1s' : -0.4606382304,
                    'C.2s' : -0.0088497203,
                    'C.2p' : -0.0719767411
                    }
        alpha = {'H.1s.H.1s' : -0.2109503575,
                 'H.1s.C.1s' : 0.1205412262,
                 'H.1s.C.2s' : -0.1985876117,
                 'H.1s.C.2p' : -0.0560772299,
                 'C.1s.C.1s' : 0.0222469851,
                 'C.1s.C.2s' : 0.0150364886,
                 'C.1s.C.2p' : 0.0001958795,
                 'C.2s.C.2s' : -0.0070553905,
                 'C.2s.C.2p' : 0.0341992648,
                 'C.2p.C.2p' : 0.0507964823,
                 'C.2p.C.2p.pi' : -0.1550251713
                 }
        beta = {'H.1s.H.1s' : 1.0219102581,
                'H.1s.C.1s' : 1.0331032280,
                'H.1s.C.2s' : 0.9580551557,
                'H.1s.C.2p' : 1.0118167832,
                'C.1s.C.1s' : 1.0001054340,
                'C.1s.C.2s' : 0.9983446940,
                'C.1s.C.2p' : 1.0005966603,
                'C.2s.C.2s' : 1.0000217409,
                'C.2s.C.2p' : 1.0062644646,
                'C.2p.C.2p' : 1.0093985074,
                'C.2p.C.2p.pi' : 0.9992768177
                }
        gamma = {'H.1s.H.1s' : 2.0000000000,
                 'H.1s.C.1s' : 2.0000000000,
                 'H.1s.C.2s' : 2.0000000000,
                 'H.1s.C.2p' : 2.0000000000,
                 'C.1s.C.1s' : 2.0000000000,
                 'C.1s.C.2s' : 2.0000000000,
                 'C.1s.C.2p' : 2.0000000000,
                 'C.2s.C.2s' : 2.0000000000,
                 'C.2s.C.2p' : 2.0000000000,
                 'C.2p.C.2p' : 2.0000000000,
                 'C.2p.C.2p.pi' : 2.0000000000
                }
        delta = {'C.2p.C.2p' : -0.161845846}
        )
"""
fc_b3lyp_hartree = """    pairwise_potential = 'u_fch_b3lyp.txt'
    fock_correction(
        diagonal = {'H.1s' : -0.5883724840,
                    'C.1s' : -2.5648517200,
                    'C.2s' : -0.4624692144,
                    'C.2p' : -0.6303377783
                    }
        alpha = {'H.1s.H.1s' : -0.5519720447,
                 'H.1s.C.1s' : -0.8317385441,
                 'H.1s.C.2s' : -0.7119280336,
                 'H.1s.C.2p' : 0.5683716947,
                 'C.1s.C.1s' : -3.8348349979,
                 'C.1s.C.2s' : -0.5758834105,
                 'C.1s.C.2p' : 2.0608200977,
                 'C.2s.C.2s' : -0.6407932567,
                 'C.2s.C.2p' : 0.4884886690,
                 'C.2p.C.2p' : 0.7918994002,
                 'C.2p.C.2p.pi' : -0.7037873600
                 }
        beta = {'H.1s.H.1s' : 0.3590609752,
                'H.1s.C.1s' : 1.0757903262,
                'H.1s.C.2s' : 0.3084011144,
                'H.1s.C.2p' : 0.5802911385,
                'C.1s.C.1s' : 2.7513970703,
                'C.1s.C.2s' : 0.4200568645,
                'C.1s.C.2p' : 1.4450215545,
                'C.2s.C.2s' : 0.1792713722,
                'C.2s.C.2p' : 0.3544647239,
                'C.2p.C.2p' : 0.7162049201,
                'C.2p.C.2p.pi' : 0.4050668194
                }
        gamma = {'H.1s.H.1s' : 1.6882461298,
                 'H.1s.C.1s' : 1.0382259985,
                 'H.1s.C.2s' : 1.6377853108,
                 'H.1s.C.2p' : 1.4559209652,
                 'C.1s.C.1s' : 1.0904536338,
                 'C.1s.C.2s' : 1.8645054364,
                 'C.1s.C.2p' : 1.0799984486,
                 'C.2s.C.2s' : 1.9416060825,
                 'C.2s.C.2p' : 1.7207351962,
                 'C.2p.C.2p' : 1.3301983525,
                 'C.2p.C.2p.pi' : 1.5630908203
                }
        delta = {'C.2p.C.2p' : 0.490157426}
        )
"""
fc_pbe_hartree = """    pairwise_potential = 'u_fch_pbe.txt'
    fock_correction(
        diagonal = {'H.1s' : -0.5323449686,
                    'C.1s' : -2.5267535898,
                    'C.2s' : -0.4625201399,
                    'C.2p' : -0.6297053844
                    }
        alpha = {'H.1s.H.1s' : -0.4976991760,
                 'H.1s.C.1s' : -0.8221174254,
                 'H.1s.C.2s' : -0.7010461862,
                 'H.1s.C.2p' : 0.5527087101,
                 'C.1s.C.1s' : -3.7752632981,
                 'C.1s.C.2s' : -0.5716227438,
                 'C.1s.C.2p' : 2.0286811779,
                 'C.2s.C.2s' : -0.6299599422,
                 'C.2s.C.2p' : 0.4813526133,
                 'C.2p.C.2p' : 0.8002777650,
                 'C.2p.C.2p.pi' : -0.6955404539
                 }
        beta = {'H.1s.H.1s' : 0.3820848496,
                'H.1s.C.1s' : 1.0979037777,
                'H.1s.C.2s' : 0.3138281709,
                'H.1s.C.2p' : 0.6262015333,
                'C.1s.C.1s' : 2.7520762987,
                'C.1s.C.2s' : 0.4222052710,
                'C.1s.C.2p' : 1.4392201545,
                'C.2s.C.2s' : 0.1802535284,
                'C.2s.C.2p' : 0.3517052748,
                'C.2p.C.2p' : 0.7154081620,
                'C.2p.C.2p.pi' : 0.4068742789
                }
        gamma = {'H.1s.H.1s' : 1.6510880288,
                 'H.1s.C.1s' : 1.0666773367,
                 'H.1s.C.2s' : 1.6553793796,
                 'H.1s.C.2p' : 1.3806133853,
                 'C.1s.C.1s' : 1.0893794405,
                 'C.1s.C.2s' : 1.8791747533,
                 'C.1s.C.2p' : 1.0853812465,
                 'C.2s.C.2s' : 1.9406329022,
                 'C.2s.C.2p' : 1.7226797406,
                 'C.2p.C.2p' : 1.3274981204,
                 'C.2p.C.2p.pi' : 1.5705040525
                }
        delta = {'C.2p.C.2p' : 0.489322087}
        )
"""
tail = """) 
"""
