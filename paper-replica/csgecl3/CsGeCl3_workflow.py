import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.spacegroup import crystal
from ase.build import bulk
from ase.io import write
from ase import Atoms


def make_struc(alat, clat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
        # Define lattice parameters
    a = alat #7.69363
    b = alat #7.69363
    c = clat #9.53490
    alpha = 90.0
    beta = 90.0
    gamma = 120.0

    # Define atom positions and types
    positions = [(0.00000, 0.00000, 0.99840),  # Cs
                 (0.00000, 0.00000, 0.51553),  # Ge
                 (0.18269, 0.36538, 0.31402)]  # Cl

    symbols = ['Cs', 'Ge', 'Cl']

    # Create ASE Atoms object
    #germanium = bulk('Ge', crystalstructure='diamond', a=5.657)
    perov = crystal(symbols=symbols, basis=positions, spacegroup=160, cellpar=[a, b, c, alpha, beta, gamma])

    #return atoms
    #perov = Atoms(symbols=symbols, scaled_positions=sc_pos, cell=lattice)
    # check how your cell looks like
    # write('s.cif', perov)
    structure = Struc(ase2struc(perov))
    return structure


def compute_energy(alat, clat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {
        "Cs": PseudoPotential(
            ptype="paw", element="Cs", functional="GGA", name="Cs.pbe-spn-kjpaw_psl.1.0.0.UPF"
        ),
        "Ge": PseudoPotential(
            ptype="paw", element="Ge", functional="GGA", name="Ge.pbe-dn-kjpaw_psl.0.2.2.UPF"
        ),
        "Cl": PseudoPotential(
            ptype="paw", element="Cl", functional="GGA", name="Cl.pbe-n-kjpaw_psl.1.0.0.UPF"
        ),
    }
    struc = make_struc(alat=alat, clat=clat)
    # fix the Pb and Ti atoms in place during relaxation
    #constraint = Constraint(atoms={"0": [0, 0, 0], "1": [0, 0, 0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=True)
    dirname = "CsGeCl3_a_{}_c_{}_ecut_{}_nk_{}".format(alat, alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ["WORKDIR"], "final-project/paper-replica/csgecl3/dft/", dirname))
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
                "ecutrho": ecut * 10,
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
    ecut = 50
    alat_list = [7.69363]
    clat_list = [9.53490]                #numpy.linspace(3.8, 4.0, 11)

    #energy_list = []
    for alat in alat_list:
        for clat in clat_list:
            print(f'Now computing DFT for a = {alat}, c = {clat} ')
            output = compute_energy(alat=alat, clat=clat, ecut=ecut, nk=nk)
            #energy_list.append(output["energy"])
            print(output)
            print('\n\n')
    
    #plt.plot(alat_list, energy_list)
    #plt.show()


if __name__ == "__main__":
    # put here the function that you actually want to run
    lattice_scan()
