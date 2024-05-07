import sys
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.io import write
from ase.spacegroup import crystal
from ase import Atoms
import os
import csv

def make_struc(alat):
    # Generate the unit cell with space group 221 (Pm-3m)
    perov = crystal(['Ba', 'Bi', 'O'], 
                        basis=[(0.5, 0.5, 0.5), (0, 0, 0), (0, 0, 0.50)],
                        spacegroup=221,
                        cellpar=[alat, alat, alat, 90, 90, 90])
    # check how your cell looks like
    # write('s.cif', perov)
    structure = Struc(ase2struc(perov))
    return structure

def compute_bands(alat, nbnd, ecut, nk, band_kpoints):
    """
    Generates a Quantum ESPRESSO input file to compute band structures.
    
    Args:
        alat (float): Lattice constant.
        nbnd (int): Number of bands.
        ecut (int): Energy cutoff.
        band_kpoints (list): List of tuples with the k-points for the band calculation.

    Returns:
        str: Path to the generated input file.
    """
    # Set up pseudopotentials as in the energy computation
    pseudopots = {
        "Ba": PseudoPotential(
            ptype="paw", element="Ba", functional="GGA", name="Ba.pbe-spn-kjpaw_psl.1.0.0.UPF"
        ),
        "Bi": PseudoPotential(
            ptype="paw", element="Bi", functional="GGA", name="Bi.pbe-dn-kjpaw_psl.1.0.0.UPF"
        ),
        "O": PseudoPotential(
            ptype="paw", element="O", functional="GGA", name="O.pbe-n-kjpaw_psl.1.0.0.UPF"
        )
    }
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=False)
    dirname = "BaBiO3_cubic_alat_{}_ecut_{}_nk_{}_bands".format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ["WORKDIR"], "final-project/final/babio3/cubic/bands", dirname))
    input_params = PWscf_inparam(
        {
            "CONTROL": {
                "calculation": "scf",
                "pseudo_dir": os.environ["QE_POTENTIALS"],
                "outdir": runpath.path,
                "disk_io": "none",
                "etot_conv_thr": 1.0e-10,
                "forc_conv_thr": 1.0e-8,
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
                "conv_thr": 1.0e-8,
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
        ncpu=2
    )
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


if __name__ == "__main__":
  # Gamma: (0 0 0)
  # X (0.5 0 0) 
  # M (0.5 0.5 0) 
  # R(0.5 0.5 0.5)
  # Back to Gamma
  # Last number is weight - can be set to zero when finding the band structure
  # Do we have to include intermediate points? Or will the code take care of that?

  band_kpoints = [(0, 0, 0, 0), (0.5, 0, 0, 0), (0.5, 0.5, 0, 0), (0.5, 0.5, 0.5, 0), (0, 0, 0, 0)]
  babio3_alat = 4.43
  ecut = 70 # Rydberg, converged by Lorenzo
  compute_bands(alat=babio3_alat, nbnd=8, ecut=ecut, nk=8, band_kpoints=band_kpoints)