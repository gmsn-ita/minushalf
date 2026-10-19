************************************************
``minushalf band-gap``
************************************************

.. code-block:: console

    $ minushalf band-gap --help                     
      
        Usage: minushalf band-gap [OPTIONS]

        Uses output files from softwares that perform ab initio calculations to
        provide the locations of VBM, CBM and the Gap value in electronvolts.The
        names of the files required for each software are listed below, it is
        worth mentioning that, for some software, their names cannot be modified.

        VASP: PROCAR, EIGENVAL, vasprun.xml QE:   {prefix}.xml (requires --input-
        name)

        Options:
        -s, --software [VASP|QE]  Specifies the software used to perform ab initio
                                calculations.  [default: VASP]
        -b, --base-path PATH      Path to folder where the relevant files are
                                located.
        -i, --indirect            Calculate indirect band gap.
        -n, --input-name PATH     Path to the input file for pw.x or projwfc.x.
                                Required for QE to determine <prefix> of output
                                files.
        --help                    Show this message and exit.

Examples
********************

To demonstrate the command usage, one calculated the positions of CBM, VBM and the Gap value of `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .



**VASP**
#############

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

**Quantum ESPRESSO (QE)**
###########################

A self consistent ``pw.x`` run was performed using the following input ``sic-2d.in`` input file.

.. code-block:: text

        &CONTROL
        calculation = 'scf'
        prefix      = 'sic_2d'
        pseudo_dir  = './pseudo/'
        outdir      = './tmp/'
        verbosity   = 'high'
        /

        &SYSTEM
        ibrav       = 0
        nat         = 2
        ntyp        = 2

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

    
After running the QE program, the :code:`minushalf band-gap` command returned the following output: 

.. code-block:: console

        $ minushalf band-gap -s QE -n sic-2d.in -i
                _                 _           _  __ 
        _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                                

        VBM: Kpoint 48, band 4 and eigenval -3.216019584303967
        CBM: Kpoint 7, band 5 and eigenval -0.6830395254853303
        Gap: 2.533eV
        _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
                        


The ``--input-name`` option (``-n``) receives the path to ``sic-2d.in``. The
command reads the ``prefix`` entry from this file to locate
the remaining output files (``<prefix>.xml``). 
The option ``-i``` was used to indicate we want the indirect gap.

References
****************

.. [1] I. Guilhon, D. S. Koda, L. G. Ferreira, M. Marques, and L. K. Teles `Phys. Rev. B 97, 045426  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045426>`_ .
.. [2] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.