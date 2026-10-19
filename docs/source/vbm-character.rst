************************************************
``minushalf vbm-character``
************************************************

.. code-block:: console

    $ minushalf vbm-character --help                     

        Usage: minushalf vbm-character [OPTIONS]

        Uses output files from softwares that perform ab initio calculations to
        discover the    last valence band (VBM) and extract, in percentage, its
        character corresponding to each    orbital type (s, p, d, ... ).  The
        names of the files required for each software are listed below, it is
        worth mentioning that their names cannot be modified.

        VASP: PROCAR, EIGENVAL, vasprun.xml 
        QE:   {prefix}.xml, {filproj}.projwfc_up (requires --input-name)

        Options:
        -s, --software [VASP|QE]  Specifies the software used to perform ab initio
                                calculations.  [default: VASP]
        -b, --base-path PATH      Path to folder where the relevant files are
                                located.
        -i, --indirect            Calculate indirect band gap.
        -n, --input-name PATH     Path to the input file for pw.x or projwfc.x.
                                Required for QE to determine <prefix> and
                                <filproj> prefixes of output files. If the pw.x
                                input file is provided,  minushalf will look for a
                                unique file with the extension *.projwfc_up in the
                                base path.
        --help                    Show this message and exit.

Examples
**************

To demonstrate the command usage, one calculated the character of the last valence band of `GaN-2d <http://2dmatpedia.org/twodmaterials/doc/2dm-2992>`_ .


**VASP**
###############

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
    
As expected for honeycomb binary materials based on III-V elements, The VBM states located at the K-point are integrally
derived from the anion :math:`p_{z}` atomic orbitals [1]_.

**Quantum ESPRESSO (QE)**
###########################

Two steps are needed to generate the needed output files: a self-consistent ``pw.x`` run using the
input file ``gan-2d.in``, followed by a projection run with ``projwfc.x`` using the
input file ``proj.in``. Below are the contents of the used input files:

.. code-block:: text

        &CONTROL
        calculation = 'scf'
        prefix      = 'gan_2d'
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
        Ga  69.723  Ga.UPF
        N   14.007  N.UPF

        CELL_PARAMETERS angstrom
        3.2180000000   0.0000000000   0.0000000000
        -1.6090000000   2.7868697494   0.0000000000
        0.0000000000   0.0000000000  20.0000000000

        ATOMIC_POSITIONS crystal
        Ga  0.3333000000  0.6666600000  0.5000000000  1  1  0
        N   0.0000000000  0.0000000000  0.5000000000  0  0  0

        K_POINTS automatic
        12 12 1 0 0 0
    
.. code-block:: text

        &PROJWFC
        outdir      = './tmp/'
        prefix      = 'gan_2d'
        filpdos     = 'mh_pdos'
        filproj     = 'mh_proj'
        lsym        = .true.
        /


After running the QE program, the :code:`minushalf vbm-character` command returned the following output: 

.. code-block:: console

        $ minushalf vbm-character -s QE -n proj.in 
                _                 _           _  __ 
        _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                                

        VBM: Kpoint 19, band 9
        |    |   d |   p |   s |
        |:---|----:|----:|----:|
        | Ga |   2 |   0 |   0 |
        | N  |   0 |  98 |   0 |
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
