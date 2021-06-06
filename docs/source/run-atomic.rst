************************************************
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
