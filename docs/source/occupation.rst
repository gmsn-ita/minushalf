************************************************
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
***************************************************

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
**************************************************************

Now, figure out a scenario where one need to generate a pseudopotential for the Ga atom with the electron medium equally divided between the orbitals :math:`p` and :math:`d`. The following command 
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

