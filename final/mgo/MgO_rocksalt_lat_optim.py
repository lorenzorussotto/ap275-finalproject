import numpy, os
import matplotlib.pyplot as plt
import sys
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.io import write, read
from ase.spacegroup import crystal
from ase.build import bulk
from ase import Atoms
import csv


def make_struc(alat):

    mgo = bulk('MgO', crystalstructure='rocksalt', a=alat)
    structure = Struc(ase2struc(mgo))
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {
        "Mg": PseudoPotential(
            ptype="paw", element="Mg", functional="GGA", name="Mg.pbe-spnl-kjpaw_psl.1.0.0.UPF"
        ),
        "O": PseudoPotential(
            ptype="paw", element="O", functional="GGA", name="O.pbe-n-kjpaw_psl.1.0.0.UPF"
        )
    }
    struc = make_struc(alat=alat)
    # fix the Pb and Ti atoms in place during relaxation
    #constraint = Constraint(atoms={"0": [0, 0, 0], "1": [0, 0, 0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=False)
    dirname = "Mg_hcp_alat_{}_ecut_{}_nk_{}".format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ["WORKDIR"], "final-project/final/mgo/rocksalt/lattice", dirname))
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
                "smearing": "mp",
                "degauss": 0.02,
            },
            "ELECTRONS": {
                "diagonalization": "david",
                "mixing_beta": 0.7,
                "conv_thr": 1e-7,
            },
            "IONS": {},
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
    alat_list = [4.15, 4.16, 4.17, 4.18, 4.19, 4.20, 4.21, 4.22, 4.23]
    
    energy_list = []

    for alat in alat_list:
        print(f'Now analyzing structure with lattice parameter = {alat}')
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        energy_list.append(round(float(output["energy"] / 2), 5))
        print(output)

    with open('MgO_lattice_optimization.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['alat', 'energy'])
        for alat, energy in zip(alat_list, energy_list):
            writer.writerow([alat, energy])


if __name__ == "__main__":
    # put here the function that you actually want to run
    lattice_scan()
