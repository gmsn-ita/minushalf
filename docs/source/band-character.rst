************************************************
``minushalf band-character``
************************************************

.. code-block:: console

    $ minushalf band-character --help                     
      
        Usage: minushalf band-character [OPTIONS] KPOINT BAND

        Uses output files from softwares that perform ab initio calculations to
        read projections in a specific kpoint band and extract, in percentage, its
        character corresponding to each orbital type (s, p, d, ... ). The   names of
        the files required for each software are listed below, it is   worth
        mentioning that, for some software, their names cannot be modified.

        VASP: PROCAR, EIGENVAL, vasprun.xml QE:   {prefix}.xml, {filproj}.projwfc_up
        (requires --input-name)

        Options:
        -s, --software [VASP|QE]  Specifies the software used to perform ab initio
                                calculations.  [default: VASP]
        -b, --base-path PATH      Path to folder where the relevant files are
                                located.
        -n, --input-name PATH     Path to the input file for pw.x or projwfc.x.
                                Required for QE to determine <prefix> and
                                <filproj> prefixes of output files. If the pw.x
                                input file is provided,  minushalf will look for a
                                unique file with the extension *.projwfc_up in the
                                base path.
        --help                    Show this message and exit.


Examples
*******************************************

To demonstrate the command usage, one calculated the character of the sixth band of the second kpoint  of `SiC-2d <http://www.2dmatpedia.org/twodmaterials/doc/2dm-2686>`_ .



**VASP**
####################

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
gradient approximation (GGA) [1]_. After running the VASP program, the :code:`minushalf band-character` command returned the following output: 

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

**Quantum ESPRESSO (QE)**
###########################

Two steps are needed to generate the needed output files: a self-consistent ``pw.x`` run using the
input file ``sic-2d.in``, followed by a projection run with ``projwfc.x`` using the
input file ``proj.in``. Below are the contents of the used input files:

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

    
.. code-block:: text

        &PROJWFC
        outdir      = './tmp/'
        prefix      = 'sic_2d'
        filpdos     = 'mh_pdos'
        filproj     = 'mh_proj'
        lsym        = .true.
        /


After running the QE program, the :code:`minushalf band-character` command returned the following output: 

.. code-block:: console

        $ minushalf band-character -s QE -n proj.in 2 6

                _                 _           _  __ 
        _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                                

        |    |   p |   s |
        |:---|----:|----:|
        | Si |  73 |   0 |
        | C  |  27 |   0 |
        _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 

The ``--input-name`` option (``-n``) receives the path to ``proj.in``. The
command reads the ``prefix`` and ``filproj`` entries from this file to locate
the remaining output files (``<prefix>.xml`` and ``<filproj>.projwfc_up``), so the 
``prefix`` entry must match the one used in the ``pw.x`` run.

In the results, band 6 of k-point 2 has 73% of Si *p* and 27% of C *p*
contribution. This differs from the VASP output because the ``KPOINT``
argument is only a position in each code's own list of k-points, so k-point 2
may not be the same point in both calculations.

References
****************

.. [1] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.