=============
Installation
=============

The minushalf CLI can be easily installed by PyPI with the following command.

.. code-block:: console

    $ pip install minushalf

Requirements
--------------
The minushalf CLI was built in order to automate the application of the DFT -1/2 [1]_ method. 
Thus, as the method requires the calculation of eigenvalues for each kpoint and band, 
it is necessary to install some software that performs ab initio calculations. 
Currently, the following softwares are supported by the program:

- VASP [2]_ [3]_
- Quantum ESPRESSO [4]_ [5]_


References
------------
.. [1] L. G. Ferreira, M. Marques, and L. K. Teles, `AIP Adv. 1, 032119 (2011) <https://doi.org/10.1063/1.3624562>`_.
.. [2] G. Kresse and J. Furthmüller, `Phys. Rev. B 54, 11169 (1996) <https://doi.org/10.1103/PhysRevB.54.11169>`_.
.. [3] G. Kresse and J. Furthmüller, `Comput. Mater. Sci. 6, 15 (1996) <https://doi.org/10.1016/0927-0256(96)00008-0>`_.
.. [4] P. Giannozzi et al., `J. Phys.: Condens. Matter 21, 395502 (2009) <https://doi.org/10.1088/0953-8984/21/39/395502>`_.
.. [5] P. Giannozzi et al., `J. Phys.: Condens. Matter 29, 465901 (2017) <https://doi.org/10.1088/1361-648X/aa8f79>`_.
