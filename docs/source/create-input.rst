************************************************
``minushalf create-input``
************************************************

This command creates the input files for the run-atomic command. Check :ref:`here <atoms_list>` the list of available atoms.

.. code-block:: console

    $ minushalf create-input --help     
    Usage: minushalf create-input [OPTIONS] CHEMICAL_SYMBOL
    
      Create the input file for the run-atomic command.

      Requires:

          CHEMICAL_SYMBOL: Chemical symbol of the atom (H, He, Na, Li...). Check the list
                           of available atoms in the docs.

      Returns:

          INP: The input file for run-atomic command

    Options:
      -e, --exchange_correlation_code [ca|wi|hl|gl|bh|pb|rp|rv|bl]
                                      Represents the functional of exchange and
                                      correlation,it can assume the following
                                      values:

                                        ca: Ceperley-Alder

                                        wi: Wigner

                                        hl: Hedin-Lundqvist

                                        gl: Gunnarson-Lundqvist

                                        bh: Von Barth-Hedin

                                        pb: PBE scheme by Perdew, Burke, and
                                        Ernzerhof
                                  
                                        rp: RPBE scheme by Hammer, Hansen, and
                                        Norskov
                                  
                                        rv: revPBE scheme by Zhang and Yang
                                  
                                        bl: BLYP (Becke-Lee-Yang-Parr) scheme
                                  
                                          [default: pb]

      -c, --calculation_code [ae]     Represents calculation code,it can assume
                                      the following values:
                                  
                                      ae: All electrons  [default: ae]

      -m, --maximum_iterations INTEGER RANGE
                                      Maximum number of iterations performed by
                                      the atomic program  [default: 100]

      -f, --filename TEXT             Name of the created file  [default: INP]
      --quiet
      --help                          Show this message and exit.
