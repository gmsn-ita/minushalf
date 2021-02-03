=========
Commands
=========

``minushalf vbm-character``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

      -b, --base-path PATH   Path to folder where the relevant files are located

      --help                 Show this message and exit.

Examples
**********

To demonstrate the command usage, one calculated the character of the last valence band of `graphene <https://materialsproject.org/materials/mp-48/>`_ .

.. container:: toggle

    .. container:: header

        **VASP**

    The following input files were used: 

    .. code-block:: xml
    
       GRAPHENE POSCAR
       1.0
       1.23 -2.130422 0.0
       1.23 2.130422  0.0
       0.00 0.00000 20
       C
       2
       Selective dynamics
       Direct
       0.0 0.0 0.5 F F F 
       0.3 0.6 0.5 T T F
       
    .. code-block:: xml

        ## INCAR
        SYSTEM = "graphene"
        ISTART = 0
        PREC = accurate
        ALGO = normal
        ISMEAR = 0; SIGMA = 0.01
        EDIFF=1.E-4
        NBANDS = 8
        NPAR = 2
        ENCUT = 350
        LORBIT = 11 ## GENERATE PROCAR
    
    .. code-block:: xml

        ## KPOINTS
        0
        Gamma
        9 9  1
        0.0 0.0 0.0

    
    The POTCAR file was used with the GGA-PBE functional of exchange and correlation. After running the VASP program,
    the :code:`minushalf vbm-character` command returned the following output: 

    .. code-block:: xml

                   _                 _           _  __ 
         _ __ ___ (_)_ __  _   _ ___| |__   __ _| |/ _|
        | '_ ` _ \| | '_ \| | | / __| '_ \ / _` | | |_ 
        | | | | | | | | | | |_| \__ \ | | | (_| | |  _|
        |_| |_| |_|_|_| |_|\__,_|___/_| |_|\__,_|_|_|  
                                               

        |    |   d |   p |   s |
        |:---|----:|----:|----:|
        | C  |   0 | 100 |   0 |

         _____ _   _ ____  
        | ____| \ | |  _ \ 
        |  _| |  \| | | | |
        | |___| |\  | |_| |
        |_____|_| \_|____/ 
    
    As expected, the character of the valence band is 100% composed
    by the carbon :math:`p` orbital.  


``minushalf cbm-character``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``minushalf band-character``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``minushalf create-input``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``minushalf run-atomic``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``minushalf occupation``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``minushalf band-gap``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``minushalf execute``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

