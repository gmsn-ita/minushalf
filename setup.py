"""
Configuration setup for python project
"""
import setuptools
from numpy.distutils.core import Extension, setup

atomic_program = Extension(
    name="minushalf.atomic_program",
    sources=[
        "minushalf/atomic_program/atm_cGuima3.f",
        "minushalf/atomic_program/atm_cGuima3.pyf"
    ],
)

setup(
    name="minushalf",
    version="1.0",
    packages=setuptools.find_packages(),
    include_package_data=True,
    author="Henrique Fernandes",
    author_email="hentt30@gmail.com",
    description=
    "CLI to provides Pre processing tools for DFT -1/2 calculations",
    license="GPL",
    install_requires=["pymatgen", "Click", "pyfiglet", "loguru"],
    entry_points="""
        [console_scripts]
        minushalf=minushalf.minushalf:minushalf
    """,
    ext_modules=[atomic_program],
)
