import numpy, os
import matplotlib.pyplot as plt
import sys
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.io import write
from ase.spacegroup import crystal
from ase import Atoms
import csv
import copy


def make_struc(alat):

    # Generate the unit cell with space group 221 (Pm-3m)
    perov = crystal(['Ba', 'Bi', 'O'], 
                        basis=[(0.5, 0.5, 0.5), (0, 0, 0), (0, 0, 0.5)],
                        spacegroup=221,
                        cellpar=[alat, alat, alat, 90, 90, 90])
    # check how your cell looks like
    # write('s.cif', perov)

    perov_f = copy.deepcopy(perov)
    perov_f[4].symbol = 'F'
    #perov_f[4].position = (0.1, perov_f[4].position[1], perov_f[4].position[2])

    structure = Struc(ase2struc(perov_f))
    return structure

def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {
        "Ba": PseudoPotential(
            ptype="paw", element="Ba", functional="GGA", name="Ba.pbe-spn-kjpaw_psl.1.0.0.UPF"
        ),
        "Bi": PseudoPotential(
            ptype="paw", element="Bi", functional="GGA", name="Bi.pbe-dn-kjpaw_psl.1.0.0.UPF"
        ),
        "O": PseudoPotential(
            ptype="paw", element="O", functional="GGA", name="O.pbe-n-kjpaw_psl.1.0.0.UPF"
        ),
        "F": PseudoPotential(
            ptype="paw", element="F", functional="GGA", name="F.pbe-n-kjpaw_psl.1.0.0.UPF"
        )
    }
    struc = make_struc(alat=alat)
    # fix the Pb and Ti atoms in place during relaxation
    #constraint = Constraint(atoms={"0": [0, 0, 0], "1": [0, 0, 0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=False)
    dirname = "BaBiO2F_cubic_alat_{}_ecut_{}_nk_{}".format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ["WORKDIR"], "final-project/final/babio2f/cubic/lattice", dirname))
    input_params = PWscf_inparam(
        {
            "CONTROL": {
                "calculation": "relax",
                "pseudo_dir": os.environ["QE_POTENTIALS"],
                "outdir": runpath.path,
                "tstress": True,
                "tprnfor": True,
                "disk_io": "none",
            },
            "SYSTEM": {
                "ecutwfc": ecut,
                "ecutrho": ecut * 8,
                "occupations": "smearing",
                "smearing": "gaussian",
                "degauss": 0.02,
            },
            "ELECTRONS": {
                "diagonalization": "david",
                "mixing_beta": 0.7,
                "conv_thr": 1e-7,
            },
            "IONS": {
                "ion_dynamics": "bfgs"
            },
            "CELL": {},
        }
    )

    output_file = run_qe_pwscf(
        runpath=runpath,
        struc=struc,
        pseudopots=pseudopots,
        params=input_params,
        kpoints=kpts,
        #constraint=constraint,
        ncpu=2
    )
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def lattice_scan():
    nk = 5
    ecut = 70
    alat_list = [4.55, 4.56, 4.57, 4.58, 4.59, 4.60, 4.61, 4.62, 4.63, 4.64, 4.65]
    
    energy_list = []

    for alat in alat_list:
        print(f'Now analyzing structure with lattice parameter = {alat}')
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        energy_list.append(round(float(output["energy"] / 5), 5))
        print(output)

    """ with open('BaBiO2F_lattice_optimization_ext2.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['alat', 'energy'])
        for alat, energy in zip(alat_list, energy_list):
            writer.writerow([alat, energy]) """


if __name__ == "__main__":
    # put here the function that you actually want to run
    lattice_scan()
