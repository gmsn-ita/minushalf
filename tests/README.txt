VASP Test Setup Instructions

To run the VASP tests, you must provide the required POTCAR files for the following elements:

Ag
C
F
S
Sb

Place these files in the directory:

test_farm/VESP/POTFILES

Each file must be named using the following convention:

POTCAR.<element>

where <element> is written in lowercase. For example:

POTCAR.ag
POTCAR.c
POTCAR.f
POTCAR.s
POTCAR.sb

Ensure that all files are correctly named and located before running the tests.