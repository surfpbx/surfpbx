import re
from collections import Counter
from typing import Union
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

from ase.units import mol
from ase.pourbaix import (
    Pourbaix,
    PourbaixDiagram,
    Species,
    parse_formula,
    RedOx,
    CONST,
    PREDEF_ENERGIES,
)


ABSOLUTE_SHE_POTENTIAL = 4.44


def convert_refs(solid: dict,
                 solvated: dict,
                 bulk_ions: dict,
                 size: int):
    """ Convert references to Surface or Species instances. """

    solid_refs = {}
    for formula, (label, energy, wf, A) in solid.items():
        surf = Surface(formula, energy, size, label, wf, A)
        solid_refs[label] = surf

    solv_refs = {}
    for formula, energy in solvated.items():
        spec = Species.from_string(formula, energy, reduce=False)
        if formula in bulk_ions:
            spec.bulk_conc = bulk_ions[formula]
        else:
            spec.bulk_conc = None
        solv_refs[formula] = spec

    return solid_refs, solv_refs


def comp2list(complist):
    return list(list(c) for c in complist)


def comp2tuple(comptup):
    return tuple(tuple(c) for c in comptup)


def draw_slice_axes(energies, X, phases, names, ptype, colors=None, min_energy=None, figsize=(12, 10)):
    from matplotlib.patches import Rectangle
    ax = plt.figure(figsize=figsize).add_subplot(111)
    ax.add_patch(Rectangle(
        (X.min(), 0),
        np.ptp(X),
        energies.max(),
        color='#fbcdc8'
    ))

    ax.plot(X, [0 for i in range(len(X))], 'k--', zorder=len(phases)+1)
    for i in range(0, energies.shape[1]):
        values = energies[:, i]
        if min_energy is not None and min(values) >= min_energy:
            continue

        phase = phases[i]
        label = names[i]

        if ptype == 'pH':
            label += f", {phase._vector[1]:.1f} e-"

        args = {
            'lw': 2.5
        }
        if colors is not None:
            args.update({'color': colors[i]})
        ax.plot(X, values, label=label, **args)

    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(energies.min(), energies.max())
    ax.set_ylabel(r'$\Delta G$ (eV)')
    plt.legend(bbox_to_anchor=(1,1), loc='upper left', fontsize=10)
    plt.subplots_adjust(
        left=0.15, right=0.75,
        top=0.9, bottom=0.14
    )
    return ax


def slice_diagram(phases, names, value, xrange,
          mode='pH', yrange=None, npoints=2,
          reference='SHE', colors=None,
          min_energy=None, figsize=(12, 10),
          savefig=None, show=False):

    X = np.linspace(*xrange, num=npoints)
    energies = np.zeros((len(X), len(phases)))

    for i, x in enumerate(X):
        for j, p in enumerate(phases):
            if mode == 'pH':
                energies[i, j] = p.evaluate(x, value)
            elif mode == 'U':
                energies[i, j] = p.evaluate(value, x)
            else:
                raise ValueError( \
                    "Unexpected mode. Choose either 'pH' or 'U'")

    ax = draw_slice_axes(energies, X, phases, names, colors=colors,
                         ptype=mode, min_energy=min_energy, figsize=figsize)

    if mode == 'pH':
        ax.set_xlabel(f'Potential vs. {reference} (V)')
    else:
        ax.set_xlabel(f'pH')

    if yrange:
        ax.set_ylim(yrange)
    if savefig:
        plt.savefig(savefig)
    if show:
        plt.show()

    return energies, ax


def format_label(label):
    return label


def format_formula_latex(formula):
    label = re.sub(r'(\S)([+-]+)', r'\1$^{\2}$', formula)
    label = re.sub(r'(\d+)', r'$_{\1}$', label)
    for symbol in ['+', '-']:
        count = label.count(symbol)
        if count > 1:
            label = label.replace(count * symbol, f'{count}{symbol}')
        if count == 1:
            label = label.replace(count * symbol, symbol)
    return label


class Surface(Species):
    def __init__(self,
                 name: str,
                 energy: float,
                 size: int = 1,
                 label: Union[str, None] = None,
                 workfunction: Union[float, None] = None,
                 area: Union[float, None] = None):

        formula, charge, aq = parse_formula(name, fmt='metal')
        if aq or charge != 0:
            raise ValueError('A surface cannot be acqueous or charged')

        super().__init__(name, formula, charge, aq, energy)

        if label is not None:
            self.name = label

        self.workfunction = workfunction
        self.area = area
        self.natoms = size


class SurfRedOx(RedOx):
    def __init__(self, combo,
                 T=298.15,
                 reference='SHE',
                 correct_conc=False,
                 size=1):

        self.species = combo
        self.T = T
        self.reference = reference
        self.size = size
        self.count = Counter()
        self.DG0 = 0
        self.correct_conc = correct_conc

        alpha = CONST * T   # 0.059 eV @ T=298.15K
        const_term = 0
        pH_term = 0
        U_term = 0
        ion = None

        for spec, coef in combo:
            if spec.aq:
                ion = spec
            self.count[spec.name] = coef
            amounts = spec.balance_electrochemistry()
            self.DG0 += coef * spec.energy
            pH_term += - coef * alpha * amounts[1]
            U_term += - coef * amounts[2]

            for name, n in zip(['H2O', 'H+', 'e-'], amounts):
                self.DG0 += coef * n * PREDEF_ENERGIES[name]
                self.count[name] += coef * n

        const_corr, U_corr, pH_corr = \
            self.get_corrections(alpha, ion=ion)

        const_term = self.DG0 + const_corr
        U_term = U_term + U_corr
        pH_term = pH_term + pH_corr

        self.DG0 /= size
        self._vector = np.array(
            [const_term, U_term, pH_term], dtype=float
        ) / size

    def get_corrections(self, alpha, ion=None, A=1, V=1):
        '''Obtain the contributions to the reaction free energy from 
             a) the surface excess-corrected ion concentration
             b) the reference electrode.
           These contribution are collected into a constant, a
           pH-dependent term and a potential-dependent term.
        '''
        A *= 1.0e16      # Electrode area converted from cm2 to Ang2
        surface = self.species[1][0]
        surf_wf = surface.workfunction
    
        ref_const_corr, ref_pH_corr = self.get_ref_correction(
            self.reference, alpha
        )
        self.DG0 += ref_const_corr
        const_corr = 0
        U_corr = 0
        pH_corr = ref_pH_corr

        if ion and self.correct_conc:
            assert surface.workfunction is not None and \
                   surface.area is not None
            charge = ion.charge
            coef = self.count[ion.name]

            # Using the user-specified bulk concentration or defaulting
            # to the concentration given by the amount of vacancies
            bulk_conc = ion.bulk_conc or \
                        1 / surface.area / mol * A / V
    
            bulk_const_term = alpha * np.log10(bulk_conc)
            met_sol_Udiff = - surf_wf + ABSOLUTE_SHE_POTENTIAL
            surf_const_term = - 0.5 * charge * met_sol_Udiff

            const_corr += coef * (bulk_const_term + surf_const_term)
            U_corr += - coef * 0.5 * charge

        return const_corr, U_corr, pH_corr

    def evaluate(self, U, pH):
        '''Evaluate Nernst equation at given pH and U'''
        K, C1, C2 = self._vector
        energy = K + np.dot((C1, C2), (U, pH))
        return energy

    def get_equilibrium_potential(self, pH=0):
        return - 1 / self._vector[1] * (self._vector[0] + pH * self._vector[2])

    def get_product_label(self):
        formatted = []
        for comp, _ in self.species:
            if self.count[comp.name] > 0:
                try:
                    formatted.append(comp.get_label())
                except AttributeError:
                    formatted.append(format_formula_latex(comp.name))
        return ', '.join(formatted)

    def copy(self):
        return self.__class__(
            self.species,
            self.T,
            self.reference,
            self.correct_conc,
            self.size
        )

    def _reinit_from_list(self, lst):
        return self.__class__(
            comp2tuple(lst),
            self.T,
            self.reference,
            self.size
        )

    def __mul__(self, n):
        comps = comp2list(self.species)
        for comp in comps:
            comp[1] *= n
        return self._reinit_from_list(comps)

    def __add__(self, other):
        comps = comp2list(self.species)
        names = [comp[0].name for comp in comps]
        for comp, coef in other.species:
            if comp.name in names:
                comps[names.index(comp.name)][1] += coef
            else:
                comps.append([comp, coef])
        return self._reinit_from_list(comps)


class SurfacePourbaix(Pourbaix):

    def __init__(self,
                 target: str,
                 solid_refs: dict,
                 aqueous_refs: dict,
                 size: int = 1,
                 T: float = 298.15,
                 correct_conc: bool = False,
                 default_conc: float = 1.0e-6,
                 bulk_ions: dict = {},
                 reference: str = 'SHE'):

        self.target = target
        self.reference = reference
        self.correct_conc = correct_conc
        self.default_conc = default_conc
        self.size = size
        self.T = T

        solid, aqueous = convert_refs(solid_refs, aqueous_refs, bulk_ions, size)
        self.material = solid.pop(target)

        self.phases, phase_matrix = self.get_surface_phases(solid, aqueous)
        self._const = phase_matrix[:, 0]
        self._var = phase_matrix[:, 1:]

    def get_surface_phases(self, products, ions):
        phases = []
        phase_matrix = []
        combos = []
    
        for product in products.values():
            diff = Counter(self.material.count)
            diff.subtract(product.count)
    
            if len(+diff) == 0:
                combos.append(((self.material, -1), (product, 1)))
                continue
    
            for elem, n in (+diff).items():
                if elem in ['O', 'H'] and abs(n) != 0:
                    combos.append(((self.material, -1), (product, 1)))
                    continue
    
                for ion in ions.values():
                    if elem in ion.count:
                        nelem = ion.count[elem]
                        combos.append((
                            (self.material, -1), (product, 1), (ion, n/nelem)
                        ))
    
        for combo in combos:
            phase = SurfRedOx(
                combo,
                self.T,
                self.reference,
                self.correct_conc,
                self.size
            )
            phases.append(phase)
            phase_matrix.append(phase._vector)
    
        return phases, np.array(phase_matrix)

    def diagram(self, U=None, pH=None):
        """ Actual generation of the diagram.

        Arguments
        ---------
        U: Iterable (1D)
            Electrostatic potential values
        pH: Iterable (1D)
            pH values

        """
        pour = np.zeros((len(U), len(pH)))
        meta = pour.copy()

        for i, u in enumerate(U):
            for j, p in enumerate(pH):
                meta[i, j], pour[i, j] = self._get_pourbaix_energy(u, p)

        # Identifying the region where the target material
        # is stable and updating the diagram accordingly
        where_stable = (meta <= 0)
        pour[where_stable] = -1

        text = []
        domains = [int(i) for i in np.unique(pour)]
        for phase_id in domains:
            if phase_id == -1:
                where = where_stable
                txt = {self.material.name: 1}
            else:
                where = (pour == phase_id)
                phase = self.phases[int(phase_id)]
                txt = phase.count
            x = np.dot(where.sum(1), U) / where.sum()
            y = np.dot(where.sum(0), pH) / where.sum()
            text.append((x, y, txt))

        return SurfacePourbaixDiagram(self, U, pH, pour, meta, text, domains)


@dataclass
class SurfacePourbaixDiagram(PourbaixDiagram):
    """ Handles diagram plotting """

    def plot(self,
             cap=1.0,
             normalize=True,
             include_text=True,
             include_water=False,
             labeltype='numbers',
             cmap="RdYlGn_r",
             filename=None,
             ax=None,
             show=False):

        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        fig = ax.get_figure()

        ax, colorbar = self._draw_diagram_axes(
            cap,
            normalize,
            include_text,
            include_water,
            labeltype, cmap,
            ax=ax)

        if normalize:
            colorbar.ax.set_ylabel(r'$\Delta G_{pbx}$ (eV/f.u.)')

        if filename is not None:
            fig.savefig(filename)

        if show:
            plt.show()
