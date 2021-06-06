************************************************
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

References
****************

.. [1] I. Guilhon, D. S. Koda, L. G. Ferreira, M. Marques, and L. K. Teles `Phys. Rev. B 97, 045426  <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.045426>`_ .
.. [2] J. P. Perdew, M. Ernzerhof, and K. Burke, `J. Chem. Phys. 105, 9982 (1996) <https://doi.org/10.1063/1.472933>`_.