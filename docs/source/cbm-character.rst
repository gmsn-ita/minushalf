************************************************
``minushalf cbm-character``
************************************************

.. code-block:: console

    $ minushalf cbm-character --help                     
      
        Usage: minushalf cbm-character [OPTIONS]

        Uses output files from softwares that perform ab initio calculations to
        discover the first conduction band (CBM) and extract, in percentage, its
        character corresponding to each orbital type (s, p, d, ... ). The   names of
        the files required for each software are listed below, it is   worth
        mentioning that their names cannot be modified.

        VASP: PROCAR, EIGENVAL, vasprun.xml 
        QE:   {prefix}.xml, {filproj}.projwfc_up (requires --input-name)

        Options:
        -s, --software [VASP|QE]  Specifies the software used to perform ab initio
                                calculations.  [default: VASP]
        -b, --base-path PATH      Path to folder where the relevant files are
                                located.
        -i, --indirect            Calculate character at indirect band gap CBM
                                k-point
        -n, --input-name PATH     Path to the input file for pw.x or projwfc.x.
                                Required for QE to determine <prefix> and
                                <filproj> prefixes of output files. If the pw.x
                                input file is provided,  minushalf will look for a
                                unique file with the extension *.projwfc_up in the
                                base path.
        --help                    Show this message and exit.


Examples
***************

To demonstrate the command usage, one calculated the character of the first conduction band of `SiC-2d <http://www.2dmatpedia.org/twodmaterials/doc/2dm-2686>`_ .



**VASP**
#######################

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


After running the QE program, the :code:`minushalf cbm-character` command returned the following output: 

.. code-block:: console

        $ minushalf cbm-character -s QE -n proj.in 

        _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                                

        CBM: Kpoint 48, band 5
        |    |   p |   s |
        |:---|----:|----:|
        | Si | 100 |   0 |
        | C  |   0 |   0 |
        _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
                              

The ``--input-name`` option (``-n``) receives the path to ``proj.in``. The
command reads the ``prefix`` and ``filproj`` entries from this file to locate
the remaining output files (``<prefix>.xml`` and ``<filproj>.projwfc_up``), so the 
``prefix`` entry must match the one used in the ``pw.x`` run.

References
****************

.. [1] I. Guilhon, D. S. Koda, L. G. Ferreira, M. Marques, and L. K. Teles `Phys. Rev. B 97, 045426  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045426>`_ .
.. [2] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.