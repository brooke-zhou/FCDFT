import os
import numpy as np

def scf_result(molecule):
# 'molecule', type=str, default=[None], nargs=1, help = "name of the molecule, eg h2"
    filename = molecule + '.out'
    if os.path.isfile(filename):
        SCF_convergence = True
        with open(filename,'r') as f:
            for lc,lines in enumerate(f):
                fields = lines.split()
                if len(fields) == 0: 
                    continue
                if fields[0] == 'Warning:' and fields[1] == 'SCF' and fields[2] == 'not' and fields[3] == 'converged!':
                    SCF_convergence = False
    else:
        print('output file does not exits.')

    # return energies
    return SCF_convergence

def fc_parse_energy(molecule,orbital):
# 'molecule', type=str, default=[None], nargs=1, help = "name of the molecule, eg h2"
# 'orbital', 'True' of 'False'. to parse orbital energies (all occupied and LUMO) or not.
    filename = molecule + '.out'
    if os.path.isfile(filename):
        if orbital == 'True':
            energy = []
            mo = []
            with open(filename,'r') as f:
                found_en = False
                found_mo = False
                for lc,lines in enumerate(f):
                    fields = lines.split()
                    if len(fields) == 0: 
                        continue
                    # if fields[0] == 'Warning:' and fields[1] == 'SCF' and fields[2] == 'not' and fields[3] == 'converged!':
                    #     print 'Warning: entos SCF not converged for ' + molecule
                    if len(fields) == 3 and fields[0] == "TOTAL" and fields[1] == "ENERGY:":
                        found_en = True
                        energy.append(float(fields[2]))
                    if (fields[0] == "Orbital" and fields[1] == "energies:") or (fields[0] == "Orbital" and fields[1] == "energies"):
                        found_mo = True
                        next(f)
                        lines = next(f)
                        while len(lines.split()) == 2:
                            mo.append(float(lines.split()[1]))
                            lines = next(f)
                if not found_en:
                    energy.append(np.nan)
                if not found_mo:
                    mo.append(np.nan)
            energy = energy + mo
            # energies = np.asarray(energy)
    
        elif orbital == 'False':
            energy = 0
            with open(filename,'r') as f:
                found_en = False
                for lc,lines in enumerate(f):
                    if lines == 'Warning: SCF not converged!':
                        print 'Warning: entos SCF not converged!'
                    fields = lines.split()
                    if len(fields) == 0: 
                        continue
                    if len(fields) == 3 and fields[0] == "TOTAL" and fields[1] == "ENERGY:":
                        found_en = True
                        energy = float(fields[2])
                if not found_en:
                    energy = np.nan
            energies = np.asarray(energy)

    else:
        print('output file does not exits.')

    # return energies
    return energy

    
        