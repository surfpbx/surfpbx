import sys
sys.path.append("..")
from surfpbx import SurfacePourbaix, slice_diagram
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


def get_kind(phase):
    ads = 0
    diss = 0
    specs = [p[0] for p in phase.species]
    print(specs)
    for p in specs:
        try:
            if 'ads' in p.name:
                ads += 1
            elif 'vac' in p.name:
                diss += 1
        except AttributeError:
            continue
    if ads == 1 and diss == 0:
        return 'ads'
    elif ads == 0 and diss == 1:
        return 'diss'
    if ads == 2 and diss == 0:
        return 'ads+ads'
    if ads == 0 and diss == 2:
        return 'diss+diss'
    if ads == 1 and diss == 1:
        return 'diss+ads'


def add_redox(p1, p2):
    spec1 = p1.count.copy()
    spec2 = p2.count.copy()
    ne1 = abs(spec1['e-'])
    ne2 = abs(spec2['e-'])
    N = max(ne1, ne2)

    new = p1 * (N / ne1) + p2 * (N / ne2)
    assert np.isclose(new.count['e-'], 0.0, rtol=0.0, atol=1e-6)
    ntarget = abs(new.species[0][1])
    return new * (1 / ntarget)


def main():
    """ Raw plot of Phosphorene aqueous stability """

    solid_refs = {
        'P16':     ('clean', -0.8264884175118539, 4.6042614015589125, 61.24877814539557),
        'H2P16':   ('ads-H-25', 1.3068808446825102, 4.407578728263564, 61.24877814539557),
        'H4P16':   ('ads-H-50', 2.7010362740593905, 4.110233998932914, 61.24877814539557),
        'H6P16':   ('ads-H-75', 3.907346848337909, 4.627379673637669, 61.24877814539557),
        'H8P16':   ('ads-H-100', 5.289880051919553, 4.4350687279084, 61.24877814539557),
        'H2O2P16': ('ads-OH-25', -4.899590841140991, 4.279446796293436, 61.24877814539557),
        'H4O4P16': ('ads-OH-50', -9.914096419038728, 4.403126095828052, 61.24877814539557),
        'H6O6P16': ('ads-OH-75', -14.922336150368068, 4.729118700645865, 61.24877814539557),
        'H8O8P16': ('ads-OH-100', -20.47331838496538, 4.760755128988246, 61.24877814539557),
        'O2P16':   ('ads-O-25', -6.247802838205814, 5.1178446960949024, 61.24877814539557),
        'O4P16':   ('ads-O-50', -10.956100323859658, 5.659557898519826, 61.24877814539557),
        'O6P16':   ('ads-O-75', -14.890331245139137, 5.940032656026439, 61.24877814539557),
        'O8P16':   ('ads-O-100', -21.997765315841434, 5.814643318849478, 61.24877814539557),
        'P15':     ('vac-P', 0.8688437048337789, 5.279013218276994, 61.24877814539557)      
    }

    aq_refs = get_solvated_refs('P')

    palette = {
        'diss': 'turquoise',
        'ads': 'darkred',
        'diss+diss': 'teal',
        'ads+ads': 'navy',
        'diss+ads': 'violet'
    }

    pbx = SurfacePourbaix(
        target='clean',         # Label of target material
        solid_refs=solid_refs,  # Solid surface phases
        aqueous_refs=aq_refs,   # Aqueous phases (ions, solvated molecules)
        size=4,                 # Number of unit cells in the modeled surface phases supercells (here: 2x2)
        correct_conc=True,      # Whether to use the surface eccess correction
        reference='Pt'          # Reference electrode
    )

    phases = []
    mixed = []
    red = []
    ox = []
    colors = []

    for p in pbx.phases:
        phases.append(p)
        if p.count['e-'] < 0:
            ox.append(p)
        elif p.count['e-'] > 0:
            red.append(p)
    
    for r in red:
        for o in ox:
            new = add_redox(r, o)
            phases.append(new)
            mixed.append(new)
    
    for p in phases:
        kind = get_kind(p)
        print(kind)
        colors.append(palette[kind])
    
    names = [p.get_product_label() for p in phases]
    slice_diagram(phases, names, 0, [0, 14], mode='U', colors=colors, show=False)
    
    plt.tight_layout()
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
