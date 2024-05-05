import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output, write_pwscf_input
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.io import write, read
from ase import Atoms


def make_struc(struct_file):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
   
    lattice = alat * numpy.identity(3)
    symbols = ["Pb", "Ti", "O", "O", "O"]
    sc_pos = [
        [0, 0, 0],
        [0.5, 0.5, 0.5],
        [0, 0.5, 0.5],
        [0.5, 0, 0.5],
        [0.5, 0.5, 0],
    ]
    perov = Atoms(symbols=symbols, scaled_positions=sc_pos, cell=lattice)
    # check how your cell looks like
    # write('s.cif', perov)
     """

    atoms = read(struct_file)
    structure = Struc(ase2struc(atoms))
    return structure


def compute_energy(alat, nk, ecut, save_path, struct_file):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {
        "Ba": PseudoPotential(
            element="Ba", name="Ba.pbe-spn-kjpaw_psl.1.0.0.UPF"
        ),
        "Bi": PseudoPotential(
            element="Bi", name="Bi.pbe-dn-kjpaw_psl.1.0.0.UPF"
        ),
        "O": PseudoPotential(
            element="O",  name="O.pbe-n-kjpaw_psl.1.0.0.UPF"
        ),
        "N": PseudoPotential(
           element="N",  name="N.pbe-n-kjpaw_psl.1.0.0.UPF"
        )
    }
    struc = make_struc(struct_file)
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=False)
    runpath = Dir(path=os.path.join(save_path, "qe-input/vc-relax"))
    input_params = PWscf_inparam(
        {
            "CONTROL": {
                "calculation": "vc-relax",
                "pseudo_dir": "../potentials/",
                "outdir": "../qe-output/",
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
            "IONS": {
                "ion_dynamics":"bfgs"
            },
            "CELL": {
                "cell_dynamics":"bfgs"
            },
        }
    )

    infile = write_pwscf_input(
        runpath=runpath, 
        params=input_params, 
        struc=struc, 
        kpoints=kpts, 
        pseudopots=pseudopots)
    
    return infile.path

def lattice_scan(save_path, struct_file):
    nk = 8
    ecut = 30
    alat = 3.9

    input_file_path = compute_energy(alat=alat, ecut=ecut, nk=nk, save_path=save_path, struct_file=struct_file) 

    print(f"Input file created at: {input_file_path}")

if __name__ == "__main__":
    # put here the function that you actually want to run
    lattice_scan(save_path='.', struct_file = 'structures/Ba8Bi8O21N3_seed2023.cif')
