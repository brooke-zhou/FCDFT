def reaction(E_atom, E_high_atom, E_mol, E_high_mol):
    RXN = 0

    # H2 --> 2H
    RXN += ( (2*E_atom['H'] - E_mol['H2'][0]) - (2*E_high_atom['H'][1] - E_high_mol['H2'][0]) )**2

    if len(E_atom) > 1: # {H,C}
        # C2H2 + C2H6 --> 2 C2H4
        RXN += ( (2*E_mol['C2H4'][0] - E_mol['C2H2'][0] - E_mol['C2H6'][0]) - \
                 (2*E_high_mol['C2H4'][0] - E_high_mol['C2H2'][0] - E_high_mol['C2H6'][0]) )**2
        # C2H4 + H2 --> C2H6
        RXN += ( (E_mol['C2H6'][0] - E_mol['C2H4'][0] - E_mol['H2'][0]) - \
                 (E_high_mol['C2H6'][0] - E_high_mol['C2H4'][0] - E_high_mol['H2'][0]) )**2
        # C2H2 + H2 --> C2H4
        RXN += ( (E_mol['C2H4'][0] - E_mol['C2H2'][0] - E_mol['H2'][0]) - \
                 (E_high_mol['C2H4'][0] - E_high_mol['C2H2'][0] - E_high_mol['H2'][0]) )**2
        # C2H2 + 2H2 --> C2H6
        RXN += ( (E_mol['C2H6'][0] - E_mol['C2H2'][0] - 2*E_mol['H2'][0]) - \
                 (E_high_mol['C2H6'][0] - E_high_mol['C2H2'][0] - 2*E_high_mol['H2'][0]) )**2
        # 3 C2H2 --> C6H6
        RXN += ( (E_mol['C6H6'][0] - 3*E_mol['C2H2'][0]) - \
                 (E_high_mol['C6H6'][0] - 3*E_high_mol['C2H2'][0]) )**2
        # C6H6 + 6 H2 --> 3 C2H6
        RXN += ( (3*E_mol['C2H6'][0] - E_mol['C6H6'][0] - 6*E_mol['H2'][0]) - \
                 (3*E_high_mol['C2H6'][0] - E_high_mol['C6H6'][0] - 6*E_high_mol['H2'][0]) )**2
        # C2H4 + 2 CH4 --> 2 C2H6
        RXN += ( (2*E_mol['C2H6'][0] - E_mol['C2H4'][0] - 2*E_mol['CH4'][0]) - \
                 (2*E_high_mol['C2H6'][0] - E_high_mol['C2H4'][0] - 2*E_high_mol['CH4'][0]) )**2
        # C2H2 + 4 CH4 --> 3 C2H6
        RXN += ( (3*E_mol['C2H6'][0] - E_mol['C2H2'][0] - 4*E_mol['CH4'][0]) - \
                 (3*E_high_mol['C2H6'][0] - E_high_mol['C2H2'][0] - 4*E_high_mol['CH4'][0]) )**2
        # C2H6 + H2 --> 2 CH4
        RXN += ( (2*E_mol['CH4'][0] - E_mol['C2H6'][0] - E_mol['H2'][0]) - \
                 (2*E_high_mol['CH4'][0] - E_high_mol['C2H6'][0] - E_high_mol['H2'][0]) )**2
        # C2H4 + 2 H2 --> 2 CH4
        RXN += ( (2*E_mol['CH4'][0] - E_mol['C2H4'][0] - 2*E_mol['H2'][0]) - \
                 (2*E_high_mol['CH4'][0] - E_high_mol['C2H4'][0] - 2*E_high_mol['H2'][0]) )**2
        # C2H2 + 3 H2 --> 2 CH4
        RXN += ( (2*E_mol['CH4'][0] - E_mol['C2H2'][0] - 3*E_mol['H2'][0]) - \
                 (2*E_high_mol['CH4'][0] - E_high_mol['C2H2'][0] - 3*E_high_mol['H2'][0]) )**2
        # C6H6 + 9 H2 --> 6 CH4
        RXN += ( (6*E_mol['CH4'][0] - E_mol['C6H6'][0] - 9*E_mol['H2'][0]) - \
                 (6*E_high_mol['CH4'][0] - E_high_mol['C6H6'][0] - 9*E_high_mol['H2'][0]) )**2

    return RXN