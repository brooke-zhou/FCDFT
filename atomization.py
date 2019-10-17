
# can only treat molecules with subscript < 100 
def mol_to_atom(molecule):
    atoms = {}

    for i, ch in enumerate(molecule):
    	if ch.isupper():
    	    if i == len(molecule) - 1 or (not molecule[i+1].islower()):
    		if atoms.has_key(ch):
    	    	    atoms[ch] = atoms[ch] + 1
    	    	else:
    	    	    atoms[ch] = 1
    	    elif molecule[i+1].islower():
    	    	element = ch + molecule[i+1]
    	    	if atoms.has_key(element):
    	    	    atoms[element] = atoms[element] + 1
    	    	else:
    	    	    atoms[element] = 1
    	elif ch.isdigit():
    	    if molecule[i-1].isdigit():
    	        continue
    	    if i == len(molecule) - 1 or molecule[i+1].isupper():
    	    	count = int(ch)
    	    elif molecule[i+1].isdigit():
    	    	count = int(ch + molecule[i+1])
    	    if molecule[i-1].isupper():
    	        atoms[molecule[i-1]] = atoms[molecule[i-1]] + count - 1
    	    elif molecule[i-1].islower:
    	        atoms[molecule[i-2:i]] = atoms[molecule[i-2:i]] + count - 1

    return atoms

# not consider E_spin for C
def atomize(molecule, E_mol, E_atom):
    count_atoms = mol_to_atom(molecule)
    E_atomized = 0
    if type(E_atom['H']) == list:
    	for element in count_atoms:
    		E_atomized = E_atomized + count_atoms[element]*E_atom[element][0]
    else:
    	for element in count_atoms:
    		E_atomized = E_atomized + count_atoms[element]*E_atom[element]

    return E_atomized - E_mol[molecule][0]
