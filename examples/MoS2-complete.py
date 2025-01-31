import sys
sys.path.append('..')
from surfpbx import SurfacePourbaix
import numpy as np
import matplotlib.pyplot as plt
from ase.phasediagram import solvated


def get_solvated_refs(name):
    """ Obtain solvated references from ase.phasediagram.solvated.

    Water and protons are removed as they are dealt with by SurfacePourbaix.
    """
    ref_dct = {}
    solv = solvated(name)
    for name, energy in solv:
        if name not in ['H+(aq)', 'H2O(aq)']:
            ref_dct[name] = energy
    return ref_dct


def main():
    # Solid Reference format: formula: (label, formation energy, workfunction, surface area)
    # Values obtained from DFT(PBE) calculations
    solid_refs = {
        'Mo9S18':   ('clean', -29.461977449325076,  5.10545594464269, 79.0200664423635),
        'Mo8S18':   ('vac-Mo', -21.759922071831994, 5.610509324350785, 79.0200664423635),
        'Mo9S17':   ('vac-S', -26.476243364598503, 5.297009797741898, 79.0200664423635),
        'HMo9OS18': ('ads-OH-11', -29.865099964424445, 4.686071046717141, 79.0200664423635),
        'Mo9O3S18': ('ads-O-33', -33.087386064846555, 5.735391315195251, 79.0200664423635),
        'Mo9O6S18': ('ads-O-66', -36.20522823756875,  6.111085459105013, 79.0200664423635),
        'Mo9O9S18': ('ads-O-100', -38.90870135334032, 6.348841416727896, 79.0200664423635),
        'H3Mo9S18': ('ads-H-33', -24.01479926757392, 4.499094152828809, 79.0200664423635),
        'Mo9OS18':  ('ads-O-11', -30.69565573093314, 5.32616311135485, 79.0200664423635),
        'HMo9S18':  ('ads-H-11', -27.502887769398654,  4.751340495330843, 79.0200664423635)
    }
    
    aq_refs = get_solvated_refs('MoS2')
    
    pbx = SurfacePourbaix(
        target='clean',                 # Label of target material
        solid_refs=solid_refs,          # Solid surface phases
        aqueous_refs=aq_refs,           # Aqueous phases (ions, solvated molecules)
        size=4,                         # Number of unit cells in the modeled surface phases supercells (here: 2x2)
        correct_conc=True,              # Whether to use the surface eccess correction
        bulk_ions={'SO4--(aq)': 1.0},   # Ionic species already present in the solution and their molar conc.
        reference='SCE'                 # Reference electrode
    )
    
    # Listing the identified phases and corresponding reactions
    for p in pbx.phases:
        print(p.equation())
    
    # Plotting the diagram
    U = np.linspace(-2, 2, 100)
    pH = np.linspace(0, 14, 100)
    diag = pbx.diagram(U, pH)
    ax = plt.figure(figsize=(12, 6)).add_subplot(111)

    diag.plot(
        show=True,
        normalize=True,
        include_text=True,
        labeltype='numbers',
        ax=ax)


if __name__ == '__main__':
    main()
