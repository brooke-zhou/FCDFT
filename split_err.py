import math
import numpy as np
import atomization
import chemical_reactions

    # E_high_mol = {'h2':[-1.161977793146, -0.37934, 0.07688]}
    # E_high_atom = {'h':-0.497431487064}
    # E_high_pec = {'pec_h2_1':-1.15907863, 'pec_h2_3':-1.16087425, 'pec_h2_4':-1.15708724, \
    # 'pec_h2_5':-1.15136746, 'pec_h2_6':-1.14432137, 'pec_h2_7':-1.13637818, \
    # 'pec_h2_8':-1.12786336, 'pec_h2_eq-':-1.16197473, 'pec_h2_eq+':-1.16197495}

    # # debug
    # E_mol = {'H2':[-1.1, -0.1, 0.1],\
    #                 'CH4': [-40.44461829, -9.878939, -0.617891, \
    #                                 -0.338772, -0.338762, -0.338759, 0.095424], ...}
    # E_atom = {'h':-0.1}
    # E_pec = {'pec_h2_1':-1.1, 'pec_h2_3':-1.1, 'pec_h2_4':-1.1, \
    # 'pec_h2_5':-1.1, 'pec_h2_6':-1.1, 'pec_h2_7':-1.1, \
    # 'pec_h2_8':-1.1, 'pec_h2_eq-':-1.1, 'pec_h2_eq+':-1.1}

def calc_delta(E_high_mol, E_high_atom, E_high_pec, E_mol, E_atom, E_pec):
    OE = 0
    PEC = 0
    AE = 0
    RXN = 0

    for molecule in E_mol:
        if molecule[0] == 'H':
            OE = OE + np.sum((np.asarray(E_mol[molecule][1:]) - np.asarray(E_high_mol[molecule])[1:])**2)
        elif molecule == 'CH4':
            OE = OE + np.sum((np.asarray(E_mol[molecule][2:7]) - np.asarray(E_high_mol[molecule])[2:])**2)
        elif molecule == 'C2H2':
            OE = OE + np.sum((np.asarray(E_mol[molecule][3:9]) - np.asarray(E_high_mol[molecule])[3:])**2)
        elif molecule == 'C2H4':
            OE = OE + np.sum((np.asarray(E_mol[molecule][3:10]) - np.asarray(E_high_mol[molecule])[3:])**2)
        elif molecule == 'C2H6':
            OE = OE + np.sum((np.asarray(E_mol[molecule][3:11]) - np.asarray(E_high_mol[molecule])[3:])**2)
        elif molecule == 'C6H6':
            OE = OE + np.sum((np.asarray(E_mol[molecule][7:23]) - np.asarray(E_high_mol[molecule])[7:])**2)

        # PEC = PEC + (E_mol[molecule][0] - E_high_mol[molecule][0])**2

        # for FCLDA, do not use this. AE calculation is only for FCH parameterization
        # for FCLDA, coefficient for AE is 0, and E_spin is obtained after all other parameters, and use AE alone as penalty function
        # if molecule[0] == 'H':
            # AE = AE + (atomization.atomize(molecule, E_mol, E_atom) - atomization.atomize(molecule, E_high_mol, E_high_atom))**2
        # elif molecule[0] == 'C': 
            # NOT FIXED!!!!!!
            # AE = AE + (atomization.atomize(molecule, E_mol, E_atom) - atomization.atomize(molecule, E_high_mol, E_high_atom))**2

    # for molecule in E_pec:
    #     PEC = PEC + (E_pec[molecule] - E_high_pec[molecule])**2 

    RXN = chemical_reactions.reaction(E_atom, E_high_atom, E_mol, E_high_mol)

    return OE,PEC,AE,RXN

def diag_rmse(E_high_mol, E_mol):
    OE = 0
    for molecule in E_mol:
        if molecule[0] == 'H':
            OE = OE + np.sum((np.asarray(E_mol[molecule][1:2]) - np.asarray(E_high_mol[molecule])[1:2])**2)
        elif molecule == 'CH4':
            OE = OE + np.sum((np.asarray(E_mol[molecule][2:6]) - np.asarray(E_high_mol[molecule])[2:6])**2)
        elif molecule == 'C2H2':
            OE = OE + np.sum((np.asarray(E_mol[molecule][3:8]) - np.asarray(E_high_mol[molecule])[3:8])**2)
        elif molecule == 'C2H4':
            OE = OE + np.sum((np.asarray(E_mol[molecule][3:9]) - np.asarray(E_high_mol[molecule])[3:9])**2)
        elif molecule == 'C2H6':
            OE = OE + np.sum((np.asarray(E_mol[molecule][3:10]) - np.asarray(E_high_mol[molecule])[3:10])**2)
        elif molecule == 'C6H6':
            OE = OE + np.sum((np.asarray(E_mol[molecule][7:22]) - np.asarray(E_high_mol[molecule])[7:22])**2)
    return math.sqrt(OE)
