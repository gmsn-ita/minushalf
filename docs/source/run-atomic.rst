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

.. note::

   **Quantum ESPRESSO.** For Quantum ESPRESSO, the atomic program described
   above (the modified ``ATOM``) is not used. The pseudopotentials are
   corrected with the native atomic program ``ld1.x``, using the option
   ``iswitch = 4`` (LDA-1/2 correction), which requires a previously generated
   pseudopotential file. See the `ld1.x documentation
   <https://www.quantum-espresso.org/Doc/INPUT_LD1.html>`_ for all the
   available input options.

   In the current version, ``ld1.x`` using ``iswitch = 4`` does **not** work with:

   - pseudopotentials whose local potential begins at zero;
   - PAW pseudopotentials.

   The relevant input parameters are:

   - ``file_pseudo``: the original pseudopotential file;
   - ``file_pseudopw``: the file where the corrected pseudopotential is written;
   - ``configts(1)`` and ``configts(2)``: the two valence configurations used
     in the correction, where the half electron is removed in one of them;
   - ``rcutv``: the cutoff distance of the 1/2 potential, mandatory for
     ``iswitch = 4``.

   For example, to remove half an electron from the C 2p orbital:

   .. code-block:: text

       &input
         title='C',
         zed=6,
         config='[He] 2s2 2p2',
         dft='PBE'
         iswitch=4
       /
       &test
         file_pseudo='C.UPF',
         file_pseudopw='C-05.upf.temp',
         configts(1)='2s2 2p2',
         configts(2)='2s2 2p1.5',
         rcutv=2.911109744123821
       /

   To "add" half an electron to the Si 3p orbital in the conduction band, the
   occupation of 1.5 goes in the first configuration instead of the second:

   .. code-block:: text

       &input
         title='Si',
         zed=14,
         config='[Ne] 3s2 3p2',
         dft='PBE'
         iswitch=4
       /
       &test
         file_pseudo='Si.UPF',
         file_pseudopw='Si-05.upf.temp',
         configts(1)='3s2 3p1.5',
         configts(2)='3s2 3p2',
         rcutv=3.304502952789202
       /