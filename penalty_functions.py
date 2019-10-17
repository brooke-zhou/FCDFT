import os,math
import numpy as np
import write_fc_input, analyze_fc_output, split_err

def debug_err(fc_param):
    # print fc_param
    result = 0
    for i in range(0,len(fc_param)):
        result = result + (fc_param[i] + 0.3) * (fc_param[i] + 0.3)
    # print result
    return result


def H_square_err(fc_param):
    command = 'rm -r FCLDA_calc'
    os.system(command)
    # write fc file (and u file)
    write_fc_input.gen_fcdft_input(['H.xyz','H2.xyz','pec_H2_1.xyz','pec_H2_3.xyz',\
        'pec_H2_4.xyz','pec_H2_5.xyz','pec_H2_6.xyz','pec_H2_7.xyz','pec_H2_8.xyz',\
        'pec_H2_eq-.xyz','pec_H2_eq+.xyz'], \
        fc_param, 'calc', [], 'LDA', 'STO-3G','1e-10')

    # do fcdft and parse energy
    os.chdir("./FCLDA_calc")
    element = ['H']
    equil_mol = ['H2']
    pec_H2 = ['pec_H2_1','pec_H2_3', 'pec_H2_4',\
                        'pec_H2_5', 'pec_H2_6','pec_H2_7',\
                        'pec_H2_8','pec_H2_eq-','pec_H2_eq+']
    molecules = element + equil_mol + pec_H2 

    E_mol = {}
    E_atom = {}
    E_pec = {}
    for molecule in molecules:
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print('Warning: entos SCF not converged for ' + molecule + '. SCF threshold will be increased.')
            os.chdir("..")
            write_fc_input.gen_fcdft_input([molecule], fc_param, 'calc', [], 'LDA', 'STO-3G','3e-10')
            os.chdir("./LDA")
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print('Warning: SCF still not converged after increasing threshold!')
        if molecule == 'H2':
            E_mol[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'True')
        elif len(molecule) == 1:
            E_atom[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')
        else:
            E_pec[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')

    # calculate square error
    # a-OE, b-PEC, c-AE, d-RXN
    # FCH, a=1,b=10,c=1,d=0
    # FCLDA, a=1,b=10,c=0,d=1
    E_high_mol = {'H2':[-1.161977793146, -0.37934, 0.07688]}
    E_high_atom = {'H':[-0.496862054327, -0.497431487064]}
    E_high_pec = {'pec_H2_1':-1.15907863, 'pec_H2_3':-1.16087425, 'pec_H2_4':-1.15708724, \
    'pec_H2_5':-1.15136746, 'pec_H2_6':-1.14432137, 'pec_H2_7':-1.13637818, \
    'pec_H2_8':-1.12786336, 'pec_H2_eq-':-1.16197473, 'pec_H2_eq+':-1.16197495}
    a = 1
    b = 10
    c = 0
    d = 1
    OE,PEC,AE,RXN = split_err.calc_delta(E_high_mol, E_high_atom, E_high_pec, E_mol, E_atom, E_pec)
    result = (a*math.sqrt(OE) + b*math.sqrt(PEC) + c*math.sqrt(AE) + d*math.sqrt(RXN)) / (a+b+c+d)
    os.chdir("..")
    # command = 'rm -r LDA' 
    # os.system(command)
    E_mol = {}
    E_atom = {}
    E_pec = {}
    return result


def H_U_square_err(fc_param):
    command = 'rm -r FCLDA_Ucalc'
    os.system(command)
    # write fc file (and u file)
    write_fc_input.gen_fcdft_input(['H.xyz','H2.xyz','pec_H2_1.xyz','pec_H2_3.xyz',\
        'pec_H2_4.xyz','pec_H2_5.xyz','pec_H2_6.xyz','pec_H2_7.xyz','pec_H2_8.xyz',\
        'pec_H2_eq-.xyz','pec_H2_eq+.xyz'], \
        fc_param, 'Ucalc', [], 'LDA', 'STO-3G','1e-10')

    # do fcdft and parse energy
    os.chdir("./FCLDA_Ucalc")
    element = ['H']
    equil_mol = ['H2']
    pec_H2 = ['pec_H2_1','pec_H2_3', 'pec_H2_4',\
                        'pec_H2_5', 'pec_H2_6','pec_H2_7',\
                        'pec_H2_8','pec_H2_eq-','pec_H2_eq+']
    molecules = element + equil_mol + pec_H2 

    E_mol = {}
    E_atom = {}
    E_pec = {}
    for molecule in molecules:
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print 'Warning: entos SCF not converged for ' + molecule + '. SCF threshold will be increased.'
            os.chdir("..")
            write_fc_input.gen_fcdft_input([molecule], fc_param, 'Ucalc', [], 'LDA', 'STO-3G','3e-10')
            os.chdir("./LDA")
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print('Warning: SCF still not converged after increasing threshold!')
        if molecule == 'H2':
            E_mol[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'True')
        elif len(molecule) == 1:
            E_atom[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')
        else:
            E_pec[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')

    # calculate square error
    # a-OE, b-PEC, c-AE, d-RXN
    # FCH, a=1,b=10,c=1,d=0
    # FCLDA, a=1,b=10,c=0,d=1
    E_high_mol = {'H2':[-1.161977793146, -0.37934, 0.07688]}
    E_high_atom = {'H':[-0.496862054327, -0.497431487064]}
    E_high_pec = {'pec_H2_1':-1.15907863, 'pec_H2_3':-1.16087425, 'pec_H2_4':-1.15708724, \
    'pec_H2_5':-1.15136746, 'pec_H2_6':-1.14432137, 'pec_H2_7':-1.13637818, \
    'pec_H2_8':-1.12786336, 'pec_H2_eq-':-1.16197473, 'pec_H2_eq+':-1.16197495}
    a = 1
    b = 10
    c = 0
    d = 1
    OE,PEC,AE,RXN = split_err.calc_delta(E_high_mol, E_high_atom, E_high_pec, E_mol, E_atom, E_pec)
    result = (a*math.sqrt(OE) + b*math.sqrt(PEC) + c*math.sqrt(AE) + d*math.sqrt(RXN)) / (a+b+c+d)
    os.chdir("..")
    E_mol = {}
    E_atom = {}
    E_pec = {}
    return result


def H_U_only_square_err(fc_param):
    command = 'rm -r FCLDA_U_only_calc'
    os.system(command)
    # write fc file (and u file)
    write_fc_input.gen_fcdft_input(['H.xyz','H2.xyz','pec_H2_1.xyz','pec_H2_3.xyz',\
        'pec_H2_4.xyz','pec_H2_5.xyz','pec_H2_6.xyz','pec_H2_7.xyz','pec_H2_8.xyz',\
        'pec_H2_eq-.xyz','pec_H2_eq+.xyz'], \
        fc_param, 'U_only_calc', [], 'LDA', 'STO-3G','1e-10')

    # do fcdft and parse energy
    os.chdir("./FCLDA_U_only_calc")
    element = ['H']
    equil_mol = ['H2']
    pec_H2 = ['pec_H2_1','pec_H2_3', 'pec_H2_4',\
                        'pec_H2_5', 'pec_H2_6','pec_H2_7',\
                        'pec_H2_8','pec_H2_eq-','pec_H2_eq+']
    molecules = element + equil_mol + pec_H2 

    E_mol = {}
    E_atom = {}
    E_pec = {}
    for molecule in molecules:
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print 'Warning: entos SCF not converged for ' + molecule + '. SCF threshold will be increased.'
            os.chdir("..")
            write_fc_input.gen_fcdft_input([molecule], fc_param, 'Ucalc', [], 'LDA', 'STO-3G','3e-10')
            os.chdir("./LDA")
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print('Warning: SCF still not converged after increasing threshold!')
        if molecule == 'H2':
            E_mol[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'True')
        elif len(molecule) == 1:
            E_atom[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')
        else:
            E_pec[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')

    # calculate square error
    # a-OE, b-PEC, c-AE, d-RXN
    # FCH, a=1,b=10,c=1,d=0
    # FCLDA, a=1,b=10,c=0,d=1
    E_high_mol = {'H2':[-1.161977793146, -0.37934, 0.07688]}
    E_high_atom = {'H':[-0.496862054327, -0.497431487064]}
    E_high_pec = {'pec_H2_1':-1.15907863, 'pec_H2_3':-1.16087425, 'pec_H2_4':-1.15708724, \
    'pec_H2_5':-1.15136746, 'pec_H2_6':-1.14432137, 'pec_H2_7':-1.13637818, \
    'pec_H2_8':-1.12786336, 'pec_H2_eq-':-1.16197473, 'pec_H2_eq+':-1.16197495}
    a = 1
    b = 10
    c = 0
    d = 1
    OE,PEC,AE,RXN = split_err.calc_delta(E_high_mol, E_high_atom, E_high_pec, E_mol, E_atom, E_pec)
    result = (a*math.sqrt(OE) + b*math.sqrt(PEC) + c*math.sqrt(AE) + d*math.sqrt(RXN)) / (a+b+c+d)
    os.chdir("..")
    E_mol = {}
    E_atom = {}
    E_pec = {}
    return result


def HC_square_err(fc_param):
    command = 'rm -r FCLDA_calc'
    os.system(command)
    # write fc file (and u file)
    # write_fc_input.gen_fcdft_input(['H.xyz','H2.xyz','pec_H2_1.xyz','pec_H2_3.xyz',\
    #     'pec_H2_4.xyz','pec_H2_5.xyz','pec_H2_6.xyz','pec_H2_7.xyz','pec_H2_8.xyz',\
    #     'pec_H2_eq-.xyz','pec_H2_eq+.xyz',\
    #     'C.xyz','CH4.xyz','pec_CH4_1.xyz','pec_CH4_2.xyz','pec_CH4_3.xyz',\
    #     'pec_CH4_5.xyz','pec_CH4_6.xyz','pec_CH4_7.xyz','pec_CH4_8.xyz',\
    #     'pec_CH4_eq-.xyz','pec_CH4_eq+.xyz',\
    #     'C2H2.xyz','pec_C2H2_1.xyz','pec_C2H2_2.xyz','pec_C2H2_3.xyz',\
    #     'pec_C2H2_5.xyz','pec_C2H2_6.xyz','pec_C2H2_7.xyz','pec_C2H2_8.xyz',\
    #     'pec_C2H2_eq-.xyz','pec_C2H2_eq+.xyz',\
    #     'C2H4.xyz','pec_C2H4_1.xyz','pec_C2H4_2.xyz','pec_C2H4_3.xyz',\
    #     'pec_C2H4_5.xyz','pec_C2H4_6.xyz','pec_C2H4_7.xyz','pec_C2H4_8.xyz',\
    #     'pec_C2H4_eq-.xyz','pec_C2H4_eq+.xyz',\
    #     'C2H6.xyz','pec_C2H6_1.xyz','pec_C2H6_2.xyz','pec_C2H6_3.xyz',\
    #     'pec_C2H6_5.xyz','pec_C2H6_6.xyz','pec_C2H6_7.xyz','pec_C2H6_8.xyz',\
    #     'pec_C2H6_eq-.xyz','pec_C2H6_eq+.xyz',\
    #     'C6H6.xyz' ], \
    #     fc_param, 'calc', [], 'LDA', 'STO-3G','1e-10')
    write_fc_input.gen_fcdft_input(['H.xyz','H2.xyz',\
        'C.xyz','CH4.xyz','C2H2.xyz','C2H4.xyz','C2H6.xyz','C6H6.xyz' ], \
        fc_param, 'calc', [], 'LDA', 'STO-3G','1e-10')

    # do fcdft and parse energy
    os.chdir("./FCLDA_calc")
    element = ['H','C']
    equil_mol = ['H2','CH4','C2H2','C2H4','C2H6','C6H6']
    pec_H2 = ['pec_H2_1','pec_H2_3', 'pec_H2_4',\
                        'pec_H2_5', 'pec_H2_6','pec_H2_7',\
                        'pec_H2_8','pec_H2_eq-','pec_H2_eq+']
    pec_CH4 = ['pec_CH4_1','pec_CH4_2','pec_CH4_3',\
                        'pec_CH4_5','pec_CH4_6','pec_CH4_7',\
                        'pec_CH4_8','pec_CH4_eq-','pec_CH4_eq+']
    pec_C2H2 = ['pec_C2H2_1','pec_C2H2_2','pec_C2H2_3',\
                            'pec_C2H2_5','pec_C2H2_6','pec_C2H2_7',\
                            'pec_C2H2_8','pec_C2H2_eq-','pec_C2H2_eq+']
    pec_C2H4 = ['pec_C2H4_1','pec_C2H4_2','pec_C2H4_3',\
                            'pec_C2H4_5','pec_C2H4_6','pec_C2H4_7',\
                            'pec_C2H4_8','pec_C2H4_eq-','pec_C2H4_eq+']
    pec_C2H6 = ['pec_C2H6_1','pec_C2H6_2','pec_C2H6_3',\
                            'pec_C2H6_5','pec_C2H6_6','pec_C2H6_7',
                            'pec_C2H6_8','pec_C2H6_eq-','pec_C2H6_eq+']
    # molecules = element + equil_mol + pec_H2 \
    #                     + pec_CH4 + pec_C2H2 + pec_C2H4 + pec_C2H6
    molecules = ['H','H2','C','CH4','C2H2','C2H4','C2H6','C6H6'] 
    
    E_mol = {}
    E_atom = {}
    E_pec = {}
    for molecule in molecules:
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print 'Warning: entos SCF not converged for ' + molecule + '. SCF threshold will be increased.'
            os.chdir("..")
            write_fc_input.gen_fcdft_input([molecule], fc_param, 'calc', [], 'LDA', 'STO-3G','3e-10')
            os.chdir("./FCLDA_calc")
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print('Warning: SCF still not converged after increasing threshold!')
        if molecule in equil_mol :
            E_mol[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'True')
        elif molecule in element :
            E_atom[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')
        else:
            E_pec[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'False')

    # calculate square error
    # a-OE, b-PEC, c-AE, d-RXN
    # FCH, a=1,b=10,c=1,d=0
    # FCLDA, a=1,b=10,c=0,d=1
    E_high_mol = {'H2':[-1.16197779, -0.379314, 0.076827], \
                  'CH4': [-40.44461829, -9.878939, -0.617891, -0.338772, -0.338762, \
                          -0.338759, 0.095424], \
                  'C2H2': [-77.214584951, -9.905086, -9.9032, -0.664979, -0.507148, \
                           -0.434242, -0.243439, -0.243439, 0.015856], \
                  'C2H4': [-78.466093895, -9.896989, -9.896357, -0.677725, -0.514872, \
                           -0.410075, -0.362738, -0.305127, -0.231719, -0.016], \
                  'C2H6': [-79.698721668, -9.885168, -9.885058, -0.672475, -0.544752, \
                           -0.377788, -0.377786, -0.317555, -0.291621, -0.291618, 0.082051], \
                  'C6H6': [-231.935100000, -9.893093,    -9.892999, -9.892985, -9.892752, \
                           -9.892738, -9.892646, -0.765174, -0.665505, -0.665488, -0.533986, \
                           -0.533975, -0.462121, -0.404791, -0.38544,    -0.365658, -0.365636, \
                           -0.31737,    -0.293006, -0.292976, -0.217961, -0.21794,    -0.026855]}
    # [E_B3LYP/6-31G*, E_PBE/6-31G*]
    E_high_atom = {'H': [-0.496862054327, -0.497431487064], 'C': [-37.707566418831, -37.823314203936]}
    E_high_pec = {'pec_H2_1':-1.15907863, 'pec_H2_3':-1.16087425, 'pec_H2_4':-1.15708724, \
    'pec_H2_5':-1.15136746, 'pec_H2_6':-1.14432137, 'pec_H2_7':-1.13637818, \
    'pec_H2_8':-1.12786336, 'pec_H2_eq-':-1.16197473, 'pec_H2_eq+':-1.16197495, \
    'pec_CH4_1':-40.42594859,'pec_CH4_2':-40.43716954,'pec_CH4_3':-40.44258,\
    'pec_CH4_5':-40.44156794,'pec_CH4_6':-40.43713569,'pec_CH4_7':-40.43102331,'pec_CH4_8':-40.4237247,\
    'pec_CH4_eq-':-40.44373365,'pec_CH4_eq+':-40.443734240,\
    'pec_C2H2_1':-77.14765155,'pec_C2H2_2':-77.18681895,'pec_C2H2_3':-77.206849,\
    'pec_C2H2_5':-77.20749054,'pec_C2H2_6':-77.19462847,'pec_C2H2_7':-77.17611297,'pec_C2H2_8':-77.15364668,\
    'pec_C2H2_eq-':-77.21250672,'pec_C2H2_eq+':-77.21251164,\
    'pec_C2H4_1':-78.41641028,'pec_C2H4_2':-78.44330054,'pec_C2H4_3':-78.45816958,\
    'pec_C2H4_5':-78.4628654,'pec_C2H4_6':-78.45671197,'pec_C2H4_7':-78.44683525,'pec_C2H4_8':-78.43429051,\
    'pec_C2H4_eq-':-78.46423952,'pec_C2H4_eq+':-78.46424281,\
    'pec_C2H6_1':-79.68066953,'pec_C2H6_2':-79.6905873,'pec_C2H6_3':-79.69568358,\
    'pec_C2H6_5':-79.69559623,'pec_C2H6_6':-79.69196915,'pec_C2H6_7':-79.686725,'pec_C2H6_8':-79.68027971,\
    'pec_C2H6_eq-':-79.69706479,'pec_C2H6_eq+':-79.69706545 }

    a = 1
    b = 10
    c = 0
    d = 1
    OE,PEC,AE,RXN = split_err.calc_delta(E_high_mol, E_high_atom, E_high_pec, E_mol, E_atom, E_pec)
    result = (a*math.sqrt(OE) + b*math.sqrt(PEC) + c*math.sqrt(AE) + d*math.sqrt(RXN)) / (a+b+c+d)
    os.chdir("..")
    # command = 'rm -r LDA' 
    # os.system(command)
    E_mol = {}
    E_atom = {}
    E_pec = {}
    return result


def HC_diag_square_err(fc_param):
    command = 'rm -r FCLDA_calc_diag'
    os.system(command)
    # write fc file (and u file)
    write_fc_input.gen_fcdft_input(\
        ['H2.xyz','CH4.xyz','C2H2.xyz','C2H4.xyz','C2H6.xyz','C6H6.xyz' ], \
        fc_param, 'calc_diag', [], 'LDA', 'STO-3G','1e-10')
    # write_fc_input.gen_fcdft_input(['H.xyz','H2.xyz',\
    #     'C.xyz','CH4.xyz','C2H2.xyz','C2H4.xyz','C2H6.xyz','C6H6.xyz' ], \
    #     fc_param, 'calc', [], 'LDA', 'STO-3G','1e-10')

    # do fcdft and parse energy
    os.chdir("./FCLDA_calc_diag")
    equil_mol = ['H2','CH4','C2H2','C2H4','C2H6','C6H6']
    molecules = equil_mol 
    E_mol = {}
    for molecule in molecules:
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print 'Warning: entos SCF not converged for ' + molecule + '. SCF threshold will be increased.'
            os.chdir("..")
            write_fc_input.gen_fcdft_input([molecule], fc_param, 'calc_diag', [], 'LDA', 'STO-3G','3e-10')
            os.chdir("./FCLDA_calc_diag")
        command = 'entos ' + molecule + '.in' + ' > ' + molecule + '.out'
        os.system(command)
        # check SCF convergence
        SCF_convergence = analyze_fc_output.scf_result(molecule)
        if SCF_convergence == False:
            print('Warning: SCF still not converged after increasing threshold!')
        if molecule in equil_mol :
            E_mol[molecule] = analyze_fc_output.fc_parse_energy(molecule, 'True')
        else:
            print('unrecognized molecule')

    # calculate square error
    # a-OE, b-PEC, c-AE, d-RXN
    # FCH, a=1,b=10,c=1,d=0
    # FCLDA, a=1,b=10,c=0,d=1
    E_high_mol = {'H2':[-1.166148283173, -0.379449, 0.040405], \
                  'CH4': [-40.464230329435, -9.867371, -0.624588, -0.345661, \
                                  -0.345657, -0.345654, 0.025820], \
                  'C2H2': [-77.25173776, -9.896738, -9.894870, -0.680828, \
                                  -0.513817, -0.445192, -0.261575, -0.261575, -0.006920], \
                  'C2H4': [-78.501314257550, -9.887511, -9.886894, -0.689026, \
                                  -0.521054, -0.419343, -0.372942, -0.311009, -0.246843, -0.034326], \
                  'C2H6': [-79.732608143195, -9.874332, -9.874212, -0.680083, -0.551026, \
                                  -0.385141, -0.385137, -0.326533, -0.297818, -0.297814, 0.019502], \
                  'C6H6': [-232.020221970006, -9.893193, -9.892609, -9.892485, -9.891970, \
                           -9.891857, -9.891245, -0.777231, -0.675074, -0.675060, -0.542377, \
                           -0.542355, -0.471169, -0.408944, -0.398185, -0.373403, -0.373394, \
                           -0.331194, -0.300419, -0.300386, -0.230969, -0.230966, -0.041878]}

    OE = split_err.diag_rmse(E_high_mol, E_mol)
    os.chdir("..")
    E_mol = {}
    return OE


def HCO_square_err(fc_param):
    # to be continued...
    return fc_param


def linear_err(fc_param):
    # example of another kind of error
    return fc_param


    


