import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.io import write
from ase.spacegroup import crystal
from ase import Atoms


def make_struc(alat):

    # Generate the unit cell with space group 221 (Pm-3m)
    perov = crystal(['Ba', 'Bi', 'O'], 
                        basis=[(0.5, 0.5, 0.5), (0, 0, 0), (0, 0, 0.5)],
                        spacegroup=221,
                        cellpar=[alat, alat, alat, 90, 90, 90])
    # check how your cell looks like
    # write('s.cif', perov)
    structure = Struc(ase2struc(perov))
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
    }
    struc = make_struc(alat=alat)
    # fix the Pb and Ti atoms in place during relaxation
    #constraint = Constraint(atoms={"0": [0, 0, 0], "1": [0, 0, 0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=False)
    dirname = "BaBiO3_cubic_alat_{}_ecut_{}_nk_{}".format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ["WORKDIR"], "Lab4/final-project/final/babio3/cubic", dirname))
    input_params = PWscf_inparam(
        {
            "CONTROL": {
                "calculation": "scf",
                "pseudo_dir": os.environ["QE_POTENTIALS"],
                "outdir": runpath.path,
                "tstress": True,
                "tprnfor": True,
                "disk_io": "none",
            },
            "SYSTEM": {
                "ecutwfc": ecut,
                "ecutrho": ecut * 8,
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
    nk = 4
    ecut = 30
    alat_list = numpy.linspace(4.30, 4.50, 21)
    
    energy_list = []
    for alat in alat_list:
        print(f'Now analyzing structure with lattice constant = {alat}')
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        energy_list.append(output["energy"])
        print(output)

    #plt.plot(alat_list, energy_list)
    #plt.show()


if __name__ == "__main__":
    # put here the function that you actually want to run
    lattice_scan()
