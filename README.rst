##############
minushalf
##############
.. image:: https://readthedocs.org/projects/minushalf/badge/
   :target: https://minushalf.readthedocs.org
   :alt: Documentation Status

.. image:: https://img.shields.io/lgtm/alerts/g/hentt30/minushalf.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/hentt30/minushalf/alerts/
   :alt: Total alerts

.. image:: https://img.shields.io/lgtm/grade/python/g/hentt30/minushalf.svg?logo=lgtm&logoWidth=18
   :target: https://lgtm.com/projects/g/hentt30/minushalf/context:python
   :alt: Language grade: Python

.. image:: https://img.shields.io/pypi/v/minushalf.svg?style=flat-square&label=PYPI%20version
   :target: https://pypi.python.org/pypi/minushalf
   :alt: Latest version released on PyPi

.. image:: https://pepy.tech/badge/minushalf
   :target: https://pepy.tech/project/minushalf
   :alt: Number of downloads
   
The DFT -1/2 method
-------------------------------

DFT-1/2, an alternative way of referring to the LDA -1/2 _ and GGA -1/2 techniques , 
is a method that method for approximate self-energy corrections within the framework of conventional Kohn-Sham DFT 
which can be used not only with the local density approximation (LDA), but also with the generalized gradient approximation (GGA).
   
The method aims to predict energy gaps results with the same precision  as the quasiparticle correction  algorithm, considered 
the state of the art for calculating energy gap of semiconductors. In addition, the computational effort of the method 
is equivalent to the standard DFT approach and and is three orders of magnitude lower than the aforementioned GW method, which allows the technique to be applied to complex systems.
   
.. figure:: docs/source/images/dft_05_demonstration.png
   :align: center
   :width: 300
   
   Fig 1. Comparison of calculated band gaps with experiment. The red square are the SCF LDA-1/2 (standard
   LDA-1/2). The crosses are standard LDA. The small gap semiconductors are metals (negative gaps), when calculated with
   LDA. LDA-1/2 corrects the situation. The band structure calculations were made with the codes VASP  and WIEN2k. 
   
   
What is minushalf?
----------------------
   
Minushalf is a command line interface (CLI) developed by the group of semiconductor materials and nanotechnology (`GMSN <http://www.gmsn.ita.br/>`_) that aims to automate 
the application of the DFT -1/2 method. The commands available in this  CLI automate both the entire process and each of its steps in order to be 
used for several purposes.

Installation
------------------
The minushalf CLI can be easily installed by PyPI with the following command.

.. code-block:: console

    $ pip install minushalf

Requirements
--------------
The minushalf CLI was built in order to automate the application of the DFT -1/2 method. 
Thus, as the method requires the calculation of eigenvalues for each kpoint and band, 
it is necessary to install some software that performs ab initio calculations. 
Currently, the following softwares are supported by the program:

- VASP 

More information
------------------------
- `Read the docs <https://minushalf.readthedocs.io/en/latest/>`_
- `Known GMSN <http://www.gmsn.ita.br/>`_
