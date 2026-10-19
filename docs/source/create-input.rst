************************************************
``minushalf create-input``
************************************************

This command creates the ``INP`` file for a specific atom, which is the half
occupied input file used by the atomic codes. Check :ref:`here <atoms_list>`
the list of available atoms. The atomic code that uses this file depends on the
software:

- **VASP:** the atomic code is the modified version of ``ATOM``, which is
  called with the ``minushalf run-atomic`` command. See :doc:`run-atomic` for
  details.
- **Quantum ESPRESSO:** the atomic code is ``ld1.x``. See
  the note at :doc:`run-atomic` for details.
.. code-block:: console

    $ minushalf create-input --help
  Usage: minushalf create-input [OPTIONS] CHEMICAL_SYMBOL

    Create the input file for the run-atomic command.

    Requires:

        CHEMICAL_SYMBOL: Chemical symbol of the atom (H, He, Na, Li...). Check
        the list                      of available atoms in the docs

    Returns:

        INP: The input file for run-atomic command

  Options:
    -s, --software [VASP|QE]        Specifies the software used to perform ab
                                    initio calculations.  [default: VASP]
    -e, --exchange_correlation_code [ca|wi|hl|gl|bh|pb|rp|rv|bl|pz]
                                    Represents the functional of exchange and
                                    correlation, it can assume the following
                                    values:
                                    
                                      ca: Ceperley-Alder (VASP only)
                                    
                                      wi: Wigner
                                    
                                      hl: Hedin-Lundqvist
                                    
                                      gl: Gunnarson-Lundqvist
                                    
                                      bh: Von Barth-Hedin (VASP only)
                                    
                                      pb: PBE scheme by Perdew, Burke, and
                                      Ernzerhof
                                    
                                      rp: RPBE scheme by Hammer, Hansen, and
                                      Norskov
                                    
                                      rv: revPBE scheme by Zhang and Yang
                                    
                                      bl: BLYP (Becke-Lee-Yang-Parr) scheme
                                    
                                      pz: Perdew-Zunger LDA (Quantum ESPRESSO
                                      only)
                                    
                                        [default: pb]
    -c, --calculation_code [ae]     Represents calculation code,it can assume
                                    the following values:
                                    
                                    ae: All electrons  [default: ae]
    -m, --maximum_iterations INTEGER RANGE
                                    Maximum number of iterations performed by
                                    the atomic program  [default: 100;
                                    0<=x<=10000]
    -p, --potential_file TEXT       Path to the pseudopotential file to be
                                    corrected
    -f, --filename TEXT             Name of the created file  [default: INP]
    -r, --cut FLOAT                 CUT parameter  [default: 0.0]
    --quiet
    --help                          Show this message and exit.
