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
    install_requires=[
        "numpy>=1.19.4",
        "pandas>=1.1.5",
        "fortranformat>=0.2.5",
        "Click>=7.1.2",
        "pyfiglet>=0.8",
        "loguru>=0.5.3",
        "tabulate>=0.8.7",
        "pyyaml>=5.3.1",
        "scipy",
    ],
    entry_points="""
        [console_scripts]
        minushalf=minushalf.minushalf:minushalf
    """,
    ext_modules=[atomic_program],
)
