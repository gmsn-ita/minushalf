
************************************************
``minushalf execute``
************************************************

This command automates the use of the DFT -1/2. It uses the Nelder-Mead algorithm [1]_ to find
the optimal values of CUT(S) and generates a text file with all the respective CUTS and the final value 
of the gap. 

.. code-block:: console

    $ minushalf execute --help
    Usage: minushalf execute [OPTIONS]

    Uses the Nelder-Mead method to find the optimal values for the CUT(S) and,
    finally, find the corrected Gap value. This command uses external software
    to perform ab initio calculations, so it must be installed in order to
    perform the command. Check the docs for a list of the softwares supported
    by the CLI.

        Requires:

            minushalf.yaml : Parameters file. Check the docs
            for a more detailed description.

            ab_initio_files: Files needed to perform the ab initio calculations.
            They must be in the same directory as the input
            file minushalf.yaml

            potential_folder: Folder with the potential files for each atom in
            the crystal. The files must be named in the following pattern
            ${POTENTIAL_FILE_NAME}.${LOWERCASE_CHEMICAL_SYMBOL}
            This does not apply to Quantum ESPRESSO.

        Returns:

            minushalf_results.dat : File that contains the optimal
            values of the cuts and the final value of the Gap.

            corrected_valence_potfiles: Potential files corrected with optimum
            valence cuts.

            corrected_conduction_potfiles: Potential files corrected with
            optimum conduction cuts.

    Options:
    --quiet
    --help   Show this message and exit.





minushalf.yaml
**************************************
:code:`minushalf.yaml` is the input file for the command :code:`execute`, each of its tags and
default values are described below.

software tag
#####################
This tag specifies the software to perform ab initio calculations. For now,
the command supports the following values for the software tag:

- VASP (Default value)
- QE (Quantum ESPRESSO tag)

.. code-block:: yaml

    software: VASP

.. code-block:: yaml

    software: QE


vasp tag
################
The vasp tag specifies the command needed to perform first principles calculations with VASP. This tag has the following fields:

- command: Command used to perform first principles calculations. (Default: ['mpirun','vasp'])

The :code:`mpirun` command is used for convenience and can be overridden depending on the local settings of the user\'s machine. The example below 
shows an use of the vasp tag in the :code:`minushalf.yaml file`:

.. code-block:: yaml

    vasp:
        command: ['mpirun','-np','6','vasp']

qe tag
################
The qe tag specifies the commands needed to perform first principles calculations with Quantum ESPRESSO. This tag has the following fields:

- pw_command: Command used to perform self consistent ``pw.x`` calculation. (Default: ['mpirun','pw.x'])
- projwfc_command: Command used to perform projection ``projwfc.x`` calculation.
- ld1_command: Command used to run ``ld1.x`` atomic code.
- input_file: Path to the ``pw.x`` input file.


.. code-block:: yaml

    qe:
        pw_command: ['mpirun', '-np', '4', '/opt/espresso/7.5/pw.x']
        projwfc_command: ['mpirun', '-np', '4', '/opt/espresso/7.5/projwfc.x']
        ld1_command: ['mpirun', '-np', '4', '/opt/espresso/7.5/ld1.x']
        input_file: ./gan-2d.in


atomic_program tag
#########################
The atomic_program tag is a set of various informations that specifies
the settings for the atomic program execution. The informations are:

- exchange_correlation_code: Functional of exchange and correlation (Default: pb)
- calculation_code: Calculation code for the atomic program (Default: ae)
- max_iterations: Maximum number of iterations performed by the atomic program (Default: 100) 

The values that the exchange_correlation_code and calculation_code tags can assume are listed below:


``exchange_correlation_code``
--------------------------------

    - ca: Ceperley-Alder (VASP only)
    - wi: Wigner
    - hl: Hedin-Lundqvist
    - gl: Gunnarson-Lundqvist
    - bh: Von Barth-Hedin (VASP only)
    - pb: PBE scheme by Perdew, Burke, and Ernzerhof
    - rp: RPBE scheme by Hammer, Hansen, and Norskov
    - rv: revPBE scheme by Zhang and Yang
    - bl: BLYP (Becke-Lee-Yang-Parr) scheme
    - pz: Perdew-Zunger LDA (Quantum ESPRESSO only)

``calculation_code``
-----------------------
    
    - ae: All electrons



Below follows an example of the atomic_program tag in the :code:`minushalf.yaml` file:

.. code-block:: yaml

    atomic_program:
        exchange_correlation_code: wi
        calculation_code: ae
        max_iterations: 200
    
correction tag
###########################
 The correction tag specifies how the DFT -1/2
  method is executed. It contains the following parameters:

- correction_code: Code that specifies the potential correction (Default: v)
- potfiles_folder: Path to folder that holds the potential files for each atom. The files must be named in the following pattern :code:`${POTENTIAL_FILE_NAME}.${LOWERCASE_CHEMICAL_SYMBOL}` (Default: minushalf_potfiles) (Does not apply to Quantum ESPRESSO)
- amplitude: Scaling factor to be used to correct the artificially generated potential. In the vast majority of cases, the amplitude value is 1.0. However,  there are some special cases where this value needs to be adjusted [5]_.  Therefore, we recommend that you do not change this value unless you know exactly what you are doing (Default: 1.0)
- valence_cut_guess: Initial Guess for the Nelder-Mead algorithm for cut in valence correction. If not provided, the default value of :math:`0.15 + 0.84d` [6]_ will be used for each optimization, where :math:`d` is the distance of the nearest neighbor in the unit cell. (Default: :math:`0.15 + 0.84d`)
- conduction_cut_guess: Initial Guess for the Nelder-Mead algorithm for cut in conduction correction. If not provided, the default value of :math:`0.15 + 0.84d` will be used for each optimization, where :math:`d` is the distance of the nearest neighbor in the unit cell. (Default: :math:`0.15 + 0.84d`)
- tolerance: Absolute tolerance for the result of the Nelder-Mead algorithm (Default: 0.01)
- fractional_valence_treshold: :ref:`Threshold  <frac_correction>` :math:`\epsilon` for fractional valence correction (Default: 10).
- fractional_conduction_treshold: :ref:`Threshold  <frac_correction>` :math:`\epsilon` for fractional conduction correction (Default: 10). 
- overwrite_vbm: In some special cases [6]_, it is necessary to consider another band as the VBM. This tag is made for these situations. It is necessary to inform the kpoint and the band number that specifies the band location. The program immediately overwrites the old projection values and uses the new values for DFT -1/2 calculations (Default: No overwrite)
- overwrite_cbm: In some special cases [6]_, it is necessary to consider another band as the CBM. This tag is made for these situations. It is necessary to inform the kpoint and the band number that specifies the band location. The program immediately overwrites the old projection values and uses the new values for DFT -1/2 calculations (Default: No overwrite)
- inplace: This tag allows you to decide whether all calculations will be done in the root folder or not. It is recommended to pass it as :code:`True` if non-self-consistent calculations are being performed for the Gap calculation, since the program only copies the input files, the output files needed for the non-self-consistent calculation will not be considered (Default: False)
- divide_character: Factor that divides the character of each atom. It is used in cases where all the bonds of an atom are with atoms of the same chemical element, as in crystals of germanium and silicon. This factor is automatically calculated by the program, however this tag will overwrite these values. So use with caution. 
- vbm_characters: This tag allows the character values of the last valence band to be provided manually. It is recommended that this tag be used with caution as it can severely impact your results.
- cbm_characters: This tag allows the character values of the first conduction band to be provided manually. It is recommended that this tag be used with caution as it can severely impact your results.

The values that the correction_code tag can assume are listed below:



``correction_code``
---------------------------

    - v: Simple valence correction
    - vf: Fractional valence correction
    - c: Simple conduction correction
    - cf: Fractional conduction correction
    - vc: Simple valence and simple conduction corrections
    - vfc: Fractional valence and simple conduction corrections
    - vcf: Simple valence and fractional conduction corrections
    - vfcf: Fractional valence and fractional conduction corrections


The example below shows an use of correction tag in the :code:`minushalf.yaml` file:

.. code-block:: yaml

        correction:
            correction_code: vf
            potfiles_folder: ../potcar
            amplitude: 3.0
            valence_cut_guess:  [["C","p",2.0],["C","s",1.5]]  # initial guesses for each orbital that contributes to the valence band
            conduction_cut_guess:  [["Si","s",1.0],["Si","p",3.5]]
            tolerance: 0.01
            fractional_valence_treshold: 15
            fractional_conduction_treshold: 23
            overwrite_vbm: [4,9] # Kpoint and band number, respectively
            overwrite_cbm: [1,3] # Kpoint and band number, respectively
            inplace: False
	    divide_character: [["C",p,1]] # divide the p character of C with one more atom
            vbm_characters: [["C","s",34],["C","s",50]] # Overwrite the characters manually
            cbm_characters: [["C","s",34],["C","s",50]]
            



Examples
***********************************************
To demonstrate the command usage, one apply the simple valence and simple conduction correction on `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .




**VASP**
#########################
  
To execute the command, the files must be provided in the following structure:

.. code-block:: console

        .
        ├── INCAR
        ├── KPOINTS
        ├── minushalf.yaml
        ├── POSCAR
        ├── POTCAR
        └── potcars
            ├── POTCAR.c
            └── POTCAR.si
    
For the input file, the following initial settings were chosen:

.. code-block:: yaml

        software: VASP
        vasp:
            command: ['mpirun','-np','4','vasp']
  
        correction:
            correction_code: vc
            potfiles_folder: ./potcars
            valence_cut_guess: [["C","p", 3.20]]
            conduction_cut_guess: [["Si","p", 3.0]]

After executing the command, one can view the result in the file :code:`minushalf_results.dat`. the file contains information
on the values obtained in the optimization of the CUT and the resulting band energy Gap (in eV). 

    .. code-block:: xml

            Valence correction cuts:
                    (C,p):3.13 a.u
            ----------------------------------------------------------------
            Conduction correction cuts:
                    (Si,p):2.77 a.u
            ----------------------------------------------------------------
            GAP: 4.37eV

         
**Quantum ESPRESSO (QE)**
###########################

To execute the command, the files must be provided in the following structure:

.. code-block:: console

        .
        ├── sic-2d.in
        ├── minushalf.yaml
        ├── C.UPF
        └── Si.UPF

Here ``sic-2d.in`` is the ``pw.x`` input file. Unlike VASP, the pseudopotential
files keep their original names in the main folder, so always keep ``pseudo_dir  = './'`` .
For the input file, the following initial settings were chosen:

.. code-block:: text

    &CONTROL
        calculation = 'scf'
        prefix      = 'sic_2d'
        pseudo_dir  = './'
        outdir      = './tmp/'
        verbosity   = 'high'
    /

    &SYSTEM
        ibrav       = 0
        nat         = 2
        ntyp        = 2
        nbnd        = 20

        ecutwfc     = 60.0
        ecutrho     = 480.0

        occupations = 'fixed'
    /

    &ELECTRONS
        conv_thr         = 1.0d-8
        startingpot      = 'atomic'
        startingwfc      = 'atomic+random'
        diagonalization  = 'david'
    /

    ATOMIC_SPECIES
    Si  28.085  Si.UPF
    C   12.011  C.UPF

    CELL_PARAMETERS angstrom
    3.100032   -0.000007    0.000001
    -1.550022    2.684696   -0.000002
    0.000006   -0.000010   20.000000

    ATOMIC_POSITIONS crystal
    Si  0.666667  0.333335  0.295447  1  1  0
    C   0.000000  0.999998  0.295392  0  1  0

    K_POINTS automatic
    12 12 1 0 0 0

And the following ``minushalf.yaml`` file was used:

.. code-block:: yaml

        software: QE

        qe:
            pw_command: ['mpirun','-np','4','pw.x']
            projwfc_command: ['mpirun','-np','4','projwfc.x']
            ld1_command: ['mpirun','-np','4','ld1.x']
            input_file: ./sic-2d.in

        atomic_program:
            exchange_correlation_code: pb
            calculation_code: ae
            max_iterations: 100

        correction:
            correction_code: vc
            potfiles_folder: ./
            valence_cut_guess: [["C","p", 3.20]]
            conduction_cut_guess: [["Si","p", 3.0]]

After executing the command, one can view the result in the file :code:`minushalf_results.dat`. The file contains information
on the values obtained in the optimization of the CUT and the resulting band energy Gap (in eV).

    .. code-block:: xml

            Valence correction cuts:
                (C,p):3.27 a.u
            ----------------------------------------------------------------
            Conduction correction cuts:
                (Si,p):3.30 a.u
            ----------------------------------------------------------------
            GAP: 4.24eV

Results
*********
For comparison purposes, the table below shows the values obtained by the method compared with
Pure GGA, functional hybrids and GW.

.. list-table:: SiC-2D band energy gap (in eV)
   :widths: 40 40 40 40 40
   :header-rows: 1

   * - GGA
     - Hybrid
     - GW
     - DFT -1/2 (VASP)
     - DFT -1/2 (QE)
   * - 2.54
     - 3.35,3.46 [2]_
     - 4.19 [3]_,4.42 [4]_
     - 4.37
     - 4.24


References
********************

.. [1]  Nelder, John A.; R. Mead (1965). A simplex method for function minimization. `Computer Journal. 7 (4): 308–313 <https://doi.org/10.1093/comjnl/7.4.308>`_.
.. [2] Y. Rao, S. Yu, and X.-M. Duan, `Phys. Chem. Chem. Phys. 19, 17250 (2017) <https://doi.org/10.1039/C7CP02616A>`_.
.. [3] H. Sahin, S. Cahangirov, M. Topsakal, E. Bekaroglu, E. Akturk, R. T. Senger, and S. Ciraci, `Phys. Rev. B 80, 155453 (2009) <https://doi.org/10.1103/PhysRevB.80.155453>`_.
.. [4] H. C. Hsueh, G. Y. Guo, and S. G. Louie, `Phys. Rev. B 84, 085404 (2011) <https://doi.org/10.1103/PhysRevB.84.085404>`_.
.. [5] C. A. Ataide, R. R. Pelá, M. Marques, L. K. Teles, J. Furthmüller, and F. Bechstedt `Phys. Rev. B 95, 045126 – Published 17 January 2017 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.045126>`_.
.. [6] L. G. Ferreira, M. Marques, and L. K. Teles, `Phys. Rev. B 78, 125116 (2008) <http://dx.doi.org/10.1103/PhysRevB.78.125116>`_.
