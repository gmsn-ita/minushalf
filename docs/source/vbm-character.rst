************************************************
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
**************

To demonstrate the command usage, one calculated the character of the last valence band of `GaN-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2992>`_ .


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
    
As expected for honeycomb binary materials based on III-V elements, The VBM states located at the Kpoint are integrally
derived from the anion :math:`p_{z}` atomic orbitals [1]_.

References
****************

.. [1] I. Guilhon, D. S. Koda, L. G. Ferreira, M. Marques, and L. K. Teles `Phys. Rev. B 97, 045426  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045426>`_ .
.. [2] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.