"""
Configuration setup for python project
"""
from os import path
import os
import setuptools
from numpy.distutils.core import Extension, setup
import sys

extra_link_args = []
os.environ["CC"] = "g++"

if sys.platform == "win32":
    extra_link_args.append("--fcompiler=gnu95")

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

atomic_program = Extension(name="minushalf.atomic_program",
                           sources=[
                               "minushalf/atomic_program/atm_cGuima3.f",
                               "minushalf/atomic_program/atm_cGuima3.pyf"
                           ],
                           extra_link_args=extra_link_args)

setup(
    name="minushalf",
    version="1.3",
    packages=setuptools.find_packages(),
    include_package_data=True,
    author="Henrique Fernandes",
    author_email="hentt30@gmail.com",
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
        "pyyaml>=5.3.1",
        "scipy>=1.5.4",
        "aenum>=3.0.0",
    ],
    entry_points="""
        [console_scripts]
        minushalf=minushalf.minushalf:minushalf
    """,
    ext_modules=[atomic_program],
)
