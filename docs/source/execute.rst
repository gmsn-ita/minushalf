
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
    perform the command. Check the docs for an list of the softwares supported
    by the CLI.

        Requires:

            minushalf.yaml : Parameters file. Check the docs
            for a more detailed description.

            ab_initio_files: Files needed to perform the ab initio
                             calculations. They must be in the same
                             directory as the input file minushalf.yaml
            
            potential_folder: Folder with the potential files for each atom in
                              the crystal. The files must be named in the following pattern 
                              ${POTENTIAL_FILE_NAME}.${LOWERCASE_CHEMICAL_SYMBOL}

        Returns:

            minushalf_results.dat : File that contains the optimal
                                    values of the cuts and the final
                                    value of the Gap.
            
            corrected_valence_potfiles: Potential files corrected with opti-mum valence cuts.

            corrected_conduction_potfiles: Potential files corrected with optimum conduction cuts.

    Options:
    --quiet
    --help   Show this message and exit.




minushalf.yaml
**************************************
:code:`minushalf.yaml` is the input file for the command :code:`execute`, each of its tags and
default values are described below.

software tag
#####################
This tag specifies the software that to perform ab initio calculations. For now,
the command supports the following values for the software tag:

- VASP (Default value)

Currently, minushalf only supports one software, but one hope to add more soon.

.. code-block:: yaml

    software: VASP

vasp tag
################
The vasp tag specifies the command needed to perform first principles calculations. This tag has the following fields:

- command: Command used to perform first principles calculations. (Default: ['mpirun','vasp'])

The :code:`mpirun` command is used for convenience and can be overridden depending on the local settings of the user's machine. The example below 
shows an use of the vasp tag in the :code:`minushalf.yaml file`:

.. code-block:: yaml

    vasp:
        command: ['mpirun','-np','6','vasp']


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
    
    - ca: Ceperley-Alder
    - wi: Wigner
    - hl: Hedin-Lundqvist
    - gl: Gunnarson-Lundqvist
    - bh: Von Barth-Hedin
    - pb: PBE scheme by Perdew, Burke, and Ernzerhof                 
    - rp: RPBE scheme by Hammer, Hansen, and Norskov
    - rv: revPBE scheme by Zhang and Yang                            
    - bl: BLYP (Becke-Lee-Yang-Parr) scheme


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

- correction_code: Code thar specifies the potential correction (Default: v)
- potfiles_folder: Path to folder that holds the potential files for each atom. The files must be named in the following pattern :code:`${POTENTIAL_FILE_NAME}.${LOWERCASE_CHEMICAL_SYMBOL}` (Default: minushalf_potfiles)
- amplitude: Scaling factor to be used to correct the artificially generated potential. In the vast majority of cases, the amplitude value is 1.0. However,  there are some special cases where this value needs to be adjusted [5]_.  Therefore, we recommend that you do not change this value unless you know exactly what you are doing (Default: 1.0)
- valence_cut_guess: Initial Guess for the Nelder-Mead algorithm for cut in valence correction. If not provided, the default value of :math:`0.15 + 0.84d` [6]_ will be used for each optimization, where :math:`d` is the distance of the nearest neighbor in the unit cell. (Default: :math:`0.15 + 0.84d`)
- conduction_cut_guess: Initial Guess for the Nelder-Mead algorithm for cut in valence correction. If not provided, the default value of :math:`0.15 + 0.84d`  will be used will be used for each optimization, where :math:`d` is the distance of the nearest neighbor in the unit cell. (Default: :math:`0.15 + 0.84d`)
- tolerance: Absolute tolerance for the result of the Nelder-Mead algorithm (Default: 0.01)
- fractionary_valence_treshold: :ref:`Treshold  <frac_correction>` :math:`\epsilon` for fractional valence correction (Default: 10). 
- overwrite_vbm: In some special cases [6]_, it is necessary to consider another band as the VBM. This tag is made for these situations. It is necessary to inform the kpoint and the band number that specifies the band location. The program immediately overwrites the old projection values and uses the new values for DFT -1/2 calculations (Default: No overwrite)
- overwrite_cbm: In some special cases [6]_, it is necessary to consider another band as the CBM. This tag is made for these situations. It is necessary to inform the kpoint and the band number that specifies the band location. The program immediately overwrites the old projection values and uses the new values for DFT -1/2 calculations (Default: No overwrite)


The values that the correction_code tag can assume are listed below:



``correction_code``
---------------------------

    - v: Simple valence correction
    - vf: Fractional valence correction
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
            valence_cut_guess: 2.0
            conduction_cut_guess: 1.0
            tolerance: 0.01
            fractionary_valence_treshold: 15
            fractionary_conduction_treshold: 23
            overwrite_vbm: [4,9] # Kpoint and band number, respectively
            overwrite_cbm: [1,3] # Kpoint and band number, respectively



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
            valence_cut_guess: 3.20
            conduction_cut_guess: 3.0

After executing the command, one can view the result in the file :code:`minushalf_results.dat`. he file contains information
on the values obtained in the optimization of the CUT and the resulting band energy Gap (in eV). 

    .. code-block:: xml

            Valence correction cuts:
                    (C):3.13 a.u
            ----------------------------------------------------------------
            Conduction correction cuts:
                    (Si):2.77 a.u
            ----------------------------------------------------------------
            GAP: 4.37eV


For comparison purposes, the table below shows the values obtained by the method compared with
Pure GGA, functional hybrids and GW.

.. list-table:: SiC-2D band energy gap (in eV)
   :widths: 40 40 40 40
   :header-rows: 1
       
   * - GGA
     - Hybrid
     - GW
     - DFT -1/2
   * - 2.54
     - 3.35,3.46 [2]_
     - 4.19 [3]_,4.42 [4]_
     - 4.37
         
  


References
********************

.. [1]  Nelder, John A.; R. Mead (1965). A simplex method for function minimization. `Computer Journal. 7 (4): 308–313 <doi:10.1093/comjnl/7.4.308>`_.
.. [2] Y. Rao, S. Yu, and X.-M. Duan, `Phys. Chem. Chem. Phys. 19, 17250 (2017) <https://doi.org/10.1039/C7CP02616A>`_.
.. [3] H. Sahin, S. Cahangirov, M. Topsakal, E. Bekaroglu, E. Akturk, R. T. Senger, and S. Ciraci, `Phys. Rev. B 80, 155453 (2009) <https://doi.org/10.1103/PhysRevB.80.155453>`_.
.. [4] H. C. Hsueh, G. Y. Guo, and S. G. Louie, `Phys. Rev. B 84, 085404 (2011) <https://doi.org/10.1103/PhysRevB.84.085404>`_.
.. [5] C. A. Ataide, R. R. Pelá, M. Marques, L. K. Teles, J. Furthmüller, and F. Bechstedt `Phys. Rev. B 95, 045126 – Published 17 January 2017 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.045126>`_.
.. [6] L. G. Ferreira, M. Marques, and L. K. Teles, `Phys. Rev. B 78, 125116 (2008) <http://dx.doi.org/10.1103/PhysRevB.78.125116>`_.
