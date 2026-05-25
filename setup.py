"""
Configuration setup for python project
"""
import os
from os import path
import sys
import setuptools
from setuptools import Extension
from setuptools.command.build_ext import build_ext


this_directory = path.abspath(path.dirname(__file__))
ATOMIC_DIR = path.join(this_directory, "minushalf", "atomic_program")

with open(path.join(this_directory, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

class f2py_Extension(Extension):
    def __init__(self, name, sourcedir):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = path.abspath(sourcedir)


class f2py_Build(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        print("Building this shit")
        os.system(
            f"cd {ext.sourcedir} && {sys.executable} -m numpy.f2py"
            f" -c atm_cGuima3.f atm_cGuima3.pyf"
            f" --include-paths={ext.sourcedir} -m atomic_program"
        )

setuptools.setup(
    name="minushalf",
    version="1.8",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        "minushalf.atomic_program": ["*.so"],   
    },
    author="Henrique Fernandes",
    author_email="dftminushalf@gmail.com",
    description=
    "CLI to provides Pre processing tools for DFT -1/2 calculations",
    long_description=long_description,
    license="GPL",
    install_requires=[
        "pandas==3.0.3",
        "fortranformat==2.0.3",
        "click==8.1.0",
        "pyfiglet==1.0.4",
        "loguru==0.7.3",
        "tabulate==0.9.0",
        "numpy==2.4.6",
        "pyyaml==6.0.1",
        "scipy==1.17.1",
        "aenum==3.1.17",
    ],
    entry_points="""
        [console_scripts]
        minushalf=minushalf.minushalf:minushalf
    """,
    ext_modules=[f2py_Extension('minushalf.atomic_program.atomic_program', ATOMIC_DIR)],
    cmdclass=dict(build_ext=f2py_Build),
)
