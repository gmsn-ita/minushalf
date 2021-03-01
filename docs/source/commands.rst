##############
Commands
##############

``minushalf vbm-character``
************************************************

.. code-block:: console

    $ minushalf vbm-character --help                     
      
      Usage: minushalf vbm-character [OPTIONS]

      Uses output files from softwares that perform ab initio calculations to
      discover the last valence band (VBM) and extract, in percentage, its
      character corresponding to each orbital type (s, p, d, ... ). The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

      VASP: PROCAR, EIGENVAL, vasprun.xml

      Options:
      -s, --software [VASP]  Specifies the software used to perform ab initio calculations.
                             [default: VASP]

      -b, --base-path PATH   Path to folder where the relevant files are located.

      --help                 Show this message and exit.

Examples
=============

To demonstrate the command usage, one calculated the character of the last valence band of `GaN-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2992>`_ .

.. container:: toggle

    .. container:: header

        **VASP**

    The following input files were used: 

    .. code-block:: xml
    
        GaN POSCAR                                 
        1.00000000000000     
        3.2180000000000004    0.0000000000000000    0.0000000000000000
        -1.6090000000000002    2.7868697493783232    0.0000000000000000
        0.0000000000000000    0.0000000000000000   20.0000000000000000
        Ga   N 
        1     1
        Selective dynamics
        Direct
        0.3333000000000013  0.6666600000000003  0.5000000000000000   T   T   F
        0.0000000000000000  0.0000000000000000  0.5000000000000000   F   F   F
       
    .. code-block:: xml
        
        PREC = Normal
        EDIFF = 0.0001
        ENCUT = 500.0
        ISMEAR= -5
        ISTART = 0
        LREAL = .FALSE.
        LORBIT=11
    
    .. code-block:: xml

       Kpoints 
       0
       Gamma
       12 12 1
       0.0 0.0 0.0

    
    Electronic properties are investigated within the DFT by applying the Perdew-Burke-Ernzerhof (PBE) functional within the general
    gradient approximation (GGA) [2]_. After running the VASP program, the :code:`minushalf vbm-character` command returned the following output: 

    .. code-block:: console

        $ minushalf vbm-character -s VASP

        _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                               

        |    |   d |   p |   s |
        |:---|----:|----:|----:|
        | Ga |  11 |   0 |   0 |
        | N  |   0 |  89 |   0 |
        _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
    
    As expected for honeycomb binary materials based on III-V elements, The VBM states located at the Kpoint are integrally
    derived from the anion :math:`p_{z}` atomic orbitals [1]_.


``minushalf cbm-character``
************************************************

.. code-block:: console

    $ minushalf cbm-character --help                     
      
      Usage: minushalf cbm-character [OPTIONS]

      Uses output files from softwares that perform ab initio calculations to
      discover the first conduction band (CBM) and extract, in percentage, its
      character corresponding to each orbital type (s, p, d, ... ). The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

      VASP: PROCAR, EIGENVAL, vasprun.xml

      Options:
      -s, --software [VASP]  Specifies the software used to perform ab initio calculations.
                             [default: VASP]

      -b, --base-path PATH   Path to folder where the relevant files are located.

      --help                 Show this message and exit.

Examples
=============

To demonstrate the command usage, one calculated the character of the first conduction band of `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .

.. container:: toggle

    .. container:: header

        **VASP**

    The following input files were used: 

    .. code-block:: xml
    
        SiC POSCAR
        1.0
        3.100032 -0.000007 0.000001
        -1.550022 2.684696 -0.000002
        0.000006 -0.000010 20.000000
        Si C
        1 1
        Selective dynamics
        direct
        0.666667 0.333335 0.295447 T T F
        0.000000 0.999998 0.295392 F T F

       
    .. code-block:: xml
        
        PREC = Normal
        EDIFF = 0.0001
        ENCUT = 500.0
        ISMEAR= -5
        ISTART = 0
        LREAL = .FALSE.
        LORBIT=11
    
    .. code-block:: xml

       Kpoints 
       0
       Gamma
       12 12 1
       0.0 0.0 0.0

    
    Electronic properties are investigated within the DFT by applying the Perdew-Burke-Ernzerhof (PBE) functional within the general
    gradient approximation (GGA) [2]_. After running the VASP program, the :code:`minushalf cbm-character` command returned the following output: 

    .. code-block:: console

        $ minushalf cbm-character -s VASP

                   _                 _           _  __ 
         _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                               

        |    |   d |   p |   s |
        |:---|----:|----:|----:|
        | Si |   0 |  85 |   0 |
        | C  |   0 |  15 |   0 |
        _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
                   
    
    As expected for honeycomb binary materials based on the IV group, The CBM states located at the Kpoint can be associated
    with the :math:`p_{z}`  orbitals of the least electronegative element [1]_.

``minushalf band-character``
************************************************

.. code-block:: console

    $ minushalf band-character --help                     
      
      Usage: minushalf band-character [OPTIONS] KPOINT BAND

      Uses output files from softwares that perform ab initio calculations to
      read projections in a specific kpoint band and extract, in percentage,
      its   character corresponding to each orbital type (s, p, d, ... ). The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

      VASP: PROCAR, EIGENVAL, vasprun.xml

      Options:
      -s, --software [VASP]  Specifies the software used to perform ab initio calculations.
                             [default: VASP]

      -b, --base-path PATH   Path to folder where the relevant files are located.

      --help                 Show this message and exit.


Examples
===========

To demonstrate the command usage, one calculated the character of the sixth band of the second kpoint  of `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .

.. container:: toggle

    .. container:: header

        **VASP**

    The following input files were used: 

    .. code-block:: xml
    
        SiC POSCAR
        1.0
        3.100032 -0.000007 0.000001
        -1.550022 2.684696 -0.000002
        0.000006 -0.000010 20.000000
        Si C
        1 1
        Selective dynamics
        direct
        0.666667 0.333335 0.295447 T T F
        0.000000 0.999998 0.295392 F T F

       
    .. code-block:: xml
        
        PREC = Normal
        EDIFF = 0.0001
        ENCUT = 500.0
        ISMEAR= -5
        ISTART = 0
        LREAL = .FALSE.
        LORBIT=11
    
    .. code-block:: xml

       Kpoints 
       0
       Gamma
       12 12 1
       0.0 0.0 0.0

    
    Electronic properties are investigated within the DFT by applying the Perdew-Burke-Ernzerhof (PBE) functional within the general
    gradient approximation (GGA) [2]_. After running the VASP program, the :code:`minushalf band-character` command returned the following output: 

    .. code-block:: console

        $ minushalf band-character 2 6 -s VASP

                  _                 _           _  __ 
        _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                               

        |    |   d |   p |   s |
        |:---|----:|----:|----:|
        | Si |   0 |   3 |   0 |
        | C  |   0 |  97 |   0 |
         _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
                   
    As one can see, band 6 of kpoint 2 has a strong character of carbon :math:`p-type` orbitals.


``minushalf band-gap``
************************************************

.. code-block:: console

    $ minushalf band-gap --help                     
      
      Usage: minushalf band-gap [OPTIONS]

      Uses output files from softwares that perform ab initio calculations to
      provide the locations of VBM, CBM and the Gap value in electronvolts.The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

      VASP: PROCAR, EIGENVAL, vasprun.xml

      Options:
      -s, --software [VASP]  Specifies the software used to perform ab initio calculations.
                             [default: VASP]

      -b, --base-path PATH   Path to folder where the relevant files are located.

      --help                 Show this message and exit.

Examples
==========

To demonstrate the command usage, one calculated the positions of CBM, VBM and the Gap value of `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .

.. container:: toggle

    .. container:: header

        **VASP**

    The following input files were used: 

    .. code-block:: xml
    
        SiC POSCAR
        1.0
        3.100032 -0.000007 0.000001
        -1.550022 2.684696 -0.000002
        0.000006 -0.000010 20.000000
        Si C
        1 1
        Selective dynamics
        direct
        0.666667 0.333335 0.295447 T T F
        0.000000 0.999998 0.295392 F T F

       
    .. code-block:: xml
        
        PREC = Normal
        EDIFF = 0.0001
        ENCUT = 500.0
        ISMEAR= -5
        ISTART = 0
        LREAL = .FALSE.
        LORBIT=11
    
    .. code-block:: xml

       Kpoints 
       0
       Gamma
       12 12 1
       0.0 0.0 0.0 

    
    Electronic properties are investigated within the DFT by applying the Perdew-Burke-Ernzerhof (PBE) functional within the general
    gradient approximation (GGA) [2]_. After running the VASP program, the :code:`minushalf band-gap` command returned the following output: 

    .. code-block:: console

        $ minushalf band-gap -s VASP

                   _                 _           _  __ 
         _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                               

        VBM: Kpoint 48, band 4 and eigenval -3.683426
        CBM: Kpoint 68, band 5 and eigenval -1.141163
        Gap: 2.542eV
         _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
                   
                   
    As expected, the Gap found is worth 2,542eV [1]_ .


``minushalf run-atomic``
************************************************

The atomic software used in this command is a modified version of the program `ATOM <https://siesta.icmab.es/SIESTA_MATERIAL/Pseudos/atom_licence.html>`_
by professor `Luiz Guimarães Ferreira <http://lattes.cnpq.br/4694847711359239>`_. The respective modifications are listed below:

-   In this version the maximum  number of interactions ('maxit') is read, just after the valence orbitals. Thus, the input files INP.pg and INP.pt must be renamed to INP.

-   Potential was generated to be added to the pseudopotential
    given by the program. The potential to be added is in the 'adiciona' file.
    The following instruction verifies that the file exists and, if it exists, is opened and read.
    
    .. code-block:: fortran

        inquire(file='adiciona',exist=lexist)
        if(lexist)  open(unit=21,file='adiciona')
    

- Creates 'VTOTAL' file with the potential related to the Schrödinger or Dirac equation.

- Creates the psfun.Guima file with the wave functions :math:`ae`, :math:`pg` and :math:`pt`.

- The pseudopotential averages are calculated for :math:`r^{2}` e :math:`r^{4}`. Electrostatic auto energy calculation is also done to valence orbitals.

.. code-block:: console

        $ minushalf run-atomic --help
        Usage: minushalf run-atomic [OPTIONS]

        Run the atomic program. The program used is a 
        modified version of ATOM by professor Luiz Guimarães Ferreira

        Requires:

            INP: The input file for the calculation.

        Returns:

            INP.ae: A copy of the input file for the calculation.

            VTOTAL.ae: Contains the atom potential.

            OUT: Contains detailed information about the run.

            AECHARGE: Contains in four columns values of r, the “up” and “down”
            parts of the total charge density, and the total core
            charge density (the charges multiplied by 4πr^2 ).

            CHARGE: is exactly identical to AECHARGE and is generated for
            backwards compatibility.

            RHO: Like CHARGE, but without the 4πr 2 factor

            AEWFNR0...AEWFNR3: All-electron valence wavefunctions as function of
            radius, for s, p, d and f valence orbitals (0,1, 2, 3, respectively — some channels might not be available).
            They include a factor of r, the s orbitals also going to zero at the
            origin.

        Options:
        --quiet
        --help   Show this message and exit.

``minushalf occupation``
************************************************

.. code-block:: console

    $ minushalf occupation --help
    Usage: minushalf occupation [OPTIONS] ORBITAL_QUANTUM_NUMBER
                            [OCCUPATION_PERCENTUAL]

    Perform fractional occupation on the atom and generate the pseudopotential
    for this occupation. The occupation can subtract any fraction of the
    electron between 0 and 0.5, half occupation is the default.

        Requires:

            ORBITAL_QUANTUM_NUMBER: A string that defines the orbital(s) in which
            the occupation will be made, it can assume four values: (0: s | 1:
            p | 2: d | 3: f). if going to modify multiple orbitals, pass a
            string with numbers separated by commas : ("0,1,2,3").

            OCCUPATION_PERCENTUAL: A string that defines percentual of half an
            electron to be used in the occupation. The default is 100%, which
            states for 0.5e. For multiple occupations in different orbitals, pass
            a string separated by commas ("100,50,40,100"). For simplicity, to
            avoid the excessive repetition of the number 100, just replace the
            number with * ("*,30,*"). If this argument is not used, the occupation
            of half electron will be made for all orbitals passed as arguments.

            INP: Input file of the run-atomic command.

        Returns:

            INP_OCC : Input file modified for fractional occupation.

            INP.ae: A copy of the input file for the calculation.

            VTOTAL_OCC: Contains the atom potential for fractional occupation.

            OUT: Contains detailed information about the run.

            AECHARGE: Contains in four columns values of r, the “up” and “down”
            parts of the total     charge density, and the total core charge
            density (the charges multiplied by 4πr^2 ).

            CHARGE: is exactly identical to AECHARGE and is generated for
            backwards compatibility.

            RHO: Like CHARGE, but without the 4πr^2 factor

            AEWFNR0...AEWFNR3: All-electron valence wavefunctions as function of
            radius, for s, p, d, and f valence orbitals (0, 1, 2, 3,
            respectively — some channels might not be available). They include
            a factor of r, the s orbitals also going to zero at the origin.

    Options:
    --quiet
    --help   Show this message and exit.

Example of occupation in only one orbital
=============================================

Suppose one need to generate a pseudopotential for the Ga atom with the occupation of half an electron in the :math:`p` orbital. The following command 
can be used for this purpose:

.. code:: console 

    $ minushalf occupation 1 100

Where the first argument represents the azimuthal quantum number for the :math:`p` orbital and the second argument represents the fraction of half an electron
that will be used in the occupation.

Initially, only the INP input file, which is shown below, needs to be provided.

.. code-block:: xml

          ae      Ga
     n=Ga c=pb
           0.0       0.0       0.0       0.0       0.0       0.0
        5    4
        4    0     2.000     0.000
        4    1     1.000     0.000
        3    2    10.000     0.000
        4    3     0.000     0.000
    100 maxit

After running the command, the following files are created 

.. code-block:: console

        .
        ├── AECHARGE
        ├── AEWFNR0
        ├── AEWFNR1
        ├── AEWFNR2 
        ├── AEWFNR3
        ├── CHARGE
        ├── fort.5
        ├── INP.ae
        ├── INP_OCC
        ├── OUT
        ├── psfun.guima
        ├── RHO
        ├── VTOTAL0
        ├── VTOTAL2
        ├── VTOTAL3
        └── VTOTAL_OCC

Where VTOTAL_OCC represents the pseudopotential for the occupation carried out and the INP_OCC file represents the
input file with the occupation of half an electron in the :math:`p` orbital, as shown below.

.. code-block:: xml

        ae      Ga
     n=Ga c=pb
           0.0       0.0       0.0       0.0       0.0       0.0
        5    4
        4    0     2.000     0.000
        4    1     0.500     0.000
        3    2    10.000     0.000
        4    3     0.000     0.000
    100 maxit

Example of occupation in multiple orbitals
================================================

Now, imagine a scenario where one need to generate a pseudopotential for the Ga atom with the electron medium equally divided between the orbitals :math:`p` and :math:`d`. The following command 
can be used for this purpose:

.. code:: console 

    $ minushalf occupation '1,2' '50,50'

Where the first argument represents the azimuthal quantum numbers for the orbitals :math:`p` and :math:`d`, while the second argument represents the fraction of half an electron
that will be used for each orbital. As the half an electron will be shared equally between the two orbitals, the fractions chosen will be :math:`50\%` for both, which corresponds
to an occupancy of a quarter of an electron for the orbitals.

Initially, only the INP input file, which is shown below, needs to be provided.

.. code-block:: xml

          ae      Ga
     n=Ga c=pb
           0.0       0.0       0.0       0.0       0.0       0.0
        5    4
        4    0     2.000     0.000
        4    1     1.000     0.000
        3    2    10.000     0.000
        4    3     0.000     0.000
    100 maxit

After executing the command, the following files are created

.. code-block:: console

        .
        ├── AECHARGE
        ├── AEWFNR0
        ├── AEWFNR1
        ├── AEWFNR2 
        ├── AEWFNR3
        ├── CHARGE
        ├── fort.5
        ├── INP.ae
        ├── INP_OCC
        ├── OUT
        ├── psfun.guima
        ├── RHO
        ├── VTOTAL0
        ├── VTOTAL2
        ├── VTOTAL3
        └── VTOTAL_OCC

Where VTOTAL_OCC represents the pseudopotential for the occupation carried out and the INP_OCC file represents the
input file with the occupation in the :math:`p` and :math:`d` orbitals, as shown below.

.. code-block:: xml

       ae      Ga
     n=Ga c=pb
           0.0       0.0       0.0       0.0       0.0       0.0
        5    4
        4    0     2.000     0.000
        4    1     0.750     0.000
        3    2     9.750     0.000
        4    3     0.000     0.000
    100 maxit


``minushalf create-input``
************************************************

This command creates the input files for the run-atomic command. Check :ref:`here <atoms_list>` the list of available atoms.

.. code-block:: console

    $ minushalf create-input --help     
    Usage: minushalf create-input [OPTIONS] CHEMICAL_SYMBOL
    
      Create the input file for the run-atomic command.

      Requires:

          CHEMICAL_SYMBOL: Chemical symbol of the atom (H, He, Na, Li...). Check the list
                           of available atoms in the docs.

      Returns:

          INP: The input file for run-atomic command

    Options:
      -e, --exchange_correlation_code [ca|wi|hl|gl|bh|pb|rp|rv|bl]
                                      Represents the functional of exchange and
                                      correlation,it can assume the following
                                      values:

                                        ca: Ceperley-Alder

                                        wi: Wigner

                                        hl: Hedin-Lundqvist

                                        gl: Gunnarson-Lundqvist

                                        bh: Von Barth-Hedin

                                        pb: PBE scheme by Perdew, Burke, and
                                        Ernzerhof
                                  
                                        rp: RPBE scheme by Hammer, Hansen, and
                                        Norskov
                                  
                                        rv: revPBE scheme by Zhang and Yang
                                  
                                        bl: BLYP (Becke-Lee-Yang-Parr) scheme
                                  
                                          [default: pb]

      -c, --calculation_code [ae]     Represents calculation code,it can assume
                                      the following values:
                                  
                                      ae: All electrons  [default: ae]

      -m, --maximum_iterations INTEGER RANGE
                                      Maximum number of iterations performed by
                                      the atomic program  [default: 100]

      -f, --filename TEXT             Name of the created file  [default: INP]
      --quiet
      --help                          Show this message and exit.


``minushalf correct-potfile``
************************************************

.. code-block:: console

    $ minushalf correct-potfile --help          
    Usage: minushalf correct-potfile [OPTIONS]

    Generate the occupied atomic potential file used for ab initio calculations.

    Requires:

        VTOTAL.ae: pseudopotential of the atom with all electrons

        VTOTAL_OCC: pseudopotential of the occupied atom

        INP_OCC: Input file for the run-atomic command of the occupied atom

        The command also needs the potential files used by the chosen software:

            VASP: POTCAR (This name can't be changed)

    Generates:

        POTFILEcut${CUT_VALUE} (If amplitude is equal to 1.0)

        POTFILEcut${CUT_VALUE}A${AMPLITUDE_VALUE} (If amplitude is different from 1.0)

    Options:
    --quiet
    -b, --base_potfile_path PATH    Path to the folder containing the potential
                                    file

    -v, --vtotal_path PATH          Path to the pseudopotential file
                                    generated by the atomic program for the atom
                                    with all electrons.  [default: VTOTAL.ae]

    -o, --vtotal_occupied_path PATH
                                    Path to the pseudopotential file
                                    generated by the atomic program for the
                                    occupied atom.  [default: VTOTAL_OCC]

    -s, --software [VASP]           Specifies the software used to make ab initio calculations.
                                    [default: VASP]

    -c, --correction [VALENCE|CONDUCTION]
                                    Indicates whether the correction should be
                                    made in the valence band or the 
                                    conduction band.  [default: VALENCE]

    -C, --cut TEXT                  distance value used to cut the potential
                                    generated artificially by fractional atomic
                                    occupation, it can be passed in two ways:
                                  
                                    unique value : float or integer. Ex: 1.0
                                  
                                    range:  begin(float|integer):pass(float|integer):end(float|integer). Ex: 1.0:0.1:2.0  
                                    [default: 2.0]

    -a, --amplitude FLOAT RANGE     Scaling factor to be used to correct the artificially generated potential.
                                    In the vast majority of cases, the amplitude value is 1.0. However, there are some
                                    special cases where this value needs to be adjusted. Therefore, we recommend that
                                    you do not change this value unless you know exactly what you are doing  [default: 1.0]

    --help                          Show this message and exit.


To consult a case where changing the amplitude value is necessary, check the reference [7]_.


``minushalf execute``
************************************************

This command automates the use of the DFT -1/2. It uses the Nelder-Mead algorithm [3]_ to find
the optimal values of CUT(S) and generates a text file with all the respective CUTS and the final value 
of the gap. To function correctly, the command requires the following files to be provided:

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
=========================
:code:`minushalf.yaml` is the input file for the command :code:`execute`, each of its tags and
default values are described below.

software tag
-----------------
This tag specifies the software that to perform ab initio calculations. For now,
the command supports the following values for the software tag:

- VASP (Default value)

Currently, minushalf only supports one software, but one hope to add more soon.

.. code-block:: yaml

    software: VASP

vasp tag
-----------------
The vasp tag is a set of various informations that specifies the settings
for executing VASP. These informations are:

- number_of_cores: The number of cores used to run VASP. (Default: 1)
- path: entry-point for the executable (Default: vasp)

Thus, the command that runs the software is :code:`mpirun -np $ {number_of_cores} $ {path}`. The example below 
shows an use of the vasp tag in the :code:`minushalf.yaml file`:

.. code-block:: yaml

    vasp:
        number_of_cores: 4
        path: vasp_bin


atomic_program tag
---------------------
The atomic_program tag is a set of various informations that specifies
the settings for the atomic program execution. The informations are:

- exchange_correlation_code: Functional of exchange and correlation (Default: pb)
- calculation_code: Calculation code for the atomic program (Default: ae)
- max_iterations: Maximum number of iterations performed by the atomic program (Default: 100) 

The values that the exchange_correlation_code and calculation_code tags can assume are listed below:

.. container:: toggle

    .. container:: header

        ``exchange_correlation_code``
    
    - ca: Ceperley-Alder
    - wi: Wigner
    - hl: Hedin-Lundqvist
    - gl: Gunnarson-Lundqvist
    - bh: Von Barth-Hedin
    - pb: PBE scheme by Perdew, Burke, and Ernzerhof                 
    - rp: RPBE scheme by Hammer, Hansen, and Norskov
    - rv: revPBE scheme by Zhang and Yang                            
    - bl: BLYP (Becke-Lee-Yang-Parr) scheme

.. container:: toggle

    .. container:: header

        ``calculation_code``
    
    - ae: All electrons



Below follows an example of the atomic_program tag in the :code:`minushalf.yaml` file:

.. code-block:: yaml

    atomic_program:
        exchange_correlation_code: wi
        calculation_code: ae
        max_iterations: 200
    
correction tag
----------------------
 The correction tag specifies how the DFT -1/2
  method is executed. It contains the following parameters:

- correction_code: Code thar specifies the potential correction (Default: v)
- potfiles_folder: Path to folder that holds the potential files for each atom. The files must be named in the following pattern :code:`${POTENTIAL_FILE_NAME}.${LOWERCASE_CHEMICAL_SYMBOL}` (Default: minushalf_potfiles)
- amplitude: caling factor to be used to correct the artificially generated potential. 
             In the vast majority of cases, the amplitude value is 1.0. However, 
             there are some special cases where this value needs to be adjusted [7]_. 
             Therefore, we recommend that you do not change this value unless you know exactly what you are doing (Default: 1.0)
- valence_cut_guess: Initial Guess for the Nelder-Mead algorithm for cut in valence correction (Default: 3.0)
- conduction_cut_guess: Initial Guess for the Nelder-Mead algorithm for cut in valence correction (Default: 2.0)
- tolerance: Minimum level of precision for the result of the Nelder-Mead algorithm (Default: 0.01)
- fractionary_valence_treshold: :ref:`Treshold  <frac_correction>` :math:`\epsilon` for fractionary valence correction (Default: 10). 
- fractionary_conduction_treshold: :ref:`Treshold  <frac_correction>` :math:`\epsilon` for fractionary conduction correction (Default: 9).

The values that the correction_code tag can assume are listed below:

.. container:: toggle

    .. container:: header

        ``correction_code``
    
    - v: Simple valence correction
    - vf: Fractionary valence correction
    - vc: Simple valence and simple conduction corrections
    - vfc: Fractionary valence and simple conduction corrections
    - vcf: Simple valence and fractionary conduction corrections
    - vfcf: Fractionary valence and fractionary conduction corrections


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


Examples
====================
To demonstrate the command usage, one apply the simple valence and simple conduction correction on `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .


.. container:: toggle

    .. container:: header

        **VASP**
  
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
            number_of_cores: 4
  
        correction:
            correction_code: vc
            potfiles_folder: ./potcars
            valence_cut_guess: 3.20
            conduction_cut_guess: 3.0

    After executing the command, one can view the result in the file :code:`minushalf_results.dat`. he file contains information
    on the values obtained in the optimization of the CUT and the resulting band energy Gap (in eV). 

        .. code-block:: xml

            Valence correction cuts:
                    (C):3.13A
            ----------------------------------------------------------------
            Conduction correction cuts:
                    (Si):2.77A
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
             - 3.35,3.46 [4]_
             - 4.19 [5]_,4.42 [6]_
             - 4.37
         
  


References
********************

.. [1] I. Guilhon, D. S. Koda, L. G. Ferreira, M. Marques, and L. K. Teles `Phys. Rev. B 97, 045426  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045426>`_ .
.. [2] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.
.. [3]  Nelder, John A.; R. Mead (1965). A simplex method for function minimization. `Computer Journal. 7 (4): 308–313 <doi:10.1093/comjnl/7.4.308>`_.
.. [4] Y. Rao, S. Yu, and X.-M. Duan, `Phys. Chem. Chem. Phys. 19, 17250 (2017) <https://doi.org/10.1039/C7CP02616A>`_.
.. [5] H. Sahin, S. Cahangirov, M. Topsakal, E. Bekaroglu, E. Akturk, R. T. Senger, and S. Ciraci, `Phys. Rev. B 80, 155453 (2009) <https://doi.org/10.1103/PhysRevB.80.155453>`_.
.. [6] H. C. Hsueh, G. Y. Guo, and S. G. Louie, `Phys. Rev. B 84, 085404 (2011) <https://doi.org/10.1103/PhysRevB.84.085404>`_.
.. [7] C. A. Ataide, R. R. Pelá, M. Marques, L. K. Teles, J. Furthmüller, and F. Bechstedt `Phys. Rev. B 95, 045126 – Published 17 January 2017 <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.95.045126>`_.
