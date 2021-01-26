=============
Introduction
=============


DFT -1/2 technique
-----------------------
DFT-1/2, an alternative way of referring to the LDA -1/2 [1]_ [2]_ and GGA -1/2 [2]_ techniques, 
is a method that performs semiconductor band-gap calculations with precision close 
to the state of the art algorithms [2]_. These technique aim to expand the half-occupation 
technique [3]_ [4]_ [5]_, formalized by Janak's theorem, to crystals using modern exchange-correlation approaches [6]_ [7]_.

The intuition of the method comes from the premise that the energy bands of a crystal are formed by the 
overlap of atomic orbitals [8]_, especially those that make up the outmost layers. This influence can 
be quantified by the projection of the wave function in a given orbital, figure 1 shows the result of these projections in VBM and CBM for ZnO. 
Thus, corrections used 
for ionization energy calculations or electronic affinity in atoms could propagate and result in a 
Gap correction for crystalline systems, since they directly influence the energy level of the orbitals [1]_.







DFT -1/2 results
----------------------

The results obtained by the application the method has the same precision [2]_ as the GW [9]_ algorithm , considered 
the state of the art for calculating the band-gap of semiconductors. In addition, the computational complexity of the method 
is equivalent to calculating the Khon-Shan gap, which allows the technique to be applied to complex systems.


References
------------------

.. [1] L. G. Ferreira, M. Marques, and L. K. Teles, `Phys. Rev. B 78, 125116 (2008) <http://dx.doi.org/10.1103/PhysRevB.78.125116>`_.

.. [2] L. G. Ferreira, M. Marques, and L. K. Teles, `AIP Adv. 1, 032119 (2011) <https://doi.org/10.1063/1.3624562>`_.

.. [3] J.C. Slater and K. H. Johnson, `Phys. Rev. B 5, 844 (1972) <http://dx.doi.org/10.1103/PhysRevB.5.844>`_.

.. [4] J.C. Slater, `Adv. Quantum Chem. 6, 1 (1972) <http://dx.doi.org/10.1016/S0065-3276(08)60541-9>`_.

.. [5] J. C. Slater and J. H. Wood, Int. J. Quant. Chem. Suppl. 4, 3 (1971).

.. [6] J. P. Perdew and A. Zunger, `Phys. Rev. B 23, 5048 (1981) <http://dx.doi.org/10.1103/PhysRevB.23.5048>`_.

.. [7] J. P. Perdew, K. Burke, and M. Ernzerhof, `Phys. Rev. Lett. 77, 3865 (1996) <http://dx.doi.org/10.1103/PhysRevLett.77.3865>`_ .

.. [8] Holgate, Sharon Ann (2009). Understanding Solid State Physics. CRC Press. pp. 177–178. ISBN 978-1-4200-1232-3.

.. [9] G. Onida, L. Reining, and A. Rubio, `Rev. Mod. Phys. 74, 601 (2002) <http://dx.doi.org/10.1103/RevModPhys.74.601>`_.
