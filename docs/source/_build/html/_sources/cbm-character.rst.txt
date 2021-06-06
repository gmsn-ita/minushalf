************************************************
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
***************

To demonstrate the command usage, one calculated the character of the first conduction band of `SiC-2d <http://www.2dmatpedia.org/2dmaterials/doc/2dm-2686>`_ .



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

References
****************

.. [1] I. Guilhon, D. S. Koda, L. G. Ferreira, M. Marques, and L. K. Teles `Phys. Rev. B 97, 045426  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045426>`_ .
.. [2] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.