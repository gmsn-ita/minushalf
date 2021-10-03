"""
Configuration setup for python project
"""
from os import path
import setuptools
from numpy.distutils.core import Extension, setup

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

atomic_program = Extension(name="minushalf.atomic_program",
                           sources=[
                               "minushalf/atomic_program/atm_cGuima3.f",
                               "minushalf/atomic_program/atm_cGuima3.pyf",
                           ],
                           extra_compile_args=['-static'])

setup(
    name="minushalf",
    version="1.5",
    packages=setuptools.find_packages(),
    include_package_data=True,
    author="Henrique Fernandes",
    author_email="dftminushalf@gmail.com",
    description=
    "CLI to provides Pre processing tools for DFT -1/2 calculations",
    long_description=long_description,
    license="GPL",
    install_requires=[
        "numpy>=1.19.4",
        "pandas>=1.1.5",
        "fortranformat>=0.2.5",
        "Click>=7.1.2",
        "pyfiglet>=0.8",
        "loguru>=0.5.3",
        "tabulate>=0.8.7",
        "pyyaml",
        "scipy>=1.5.4",
        "aenum>=3.0.0",
    ],
    entry_points="""
        [console_scripts]
        minushalf=minushalf.minushalf:minushalf
    """,
    ext_modules=[atomic_program],
)
