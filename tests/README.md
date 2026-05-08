# VASP Test Setup Instructions

To run the VASP tests, the required POTCAR files for the following elements must be provided:

- `Ag`
- `C`
- `F`
- `S`
- `Sb`

Place these files in the following directory:

```
test_farm/VASP/POTFILES
```

Each file must follow this naming scheme:

```
POTCAR.<element>
```

where `<element>` is written in lowercase. For example:

```
POTCAR.ag
POTCAR.c
POTCAR.f
POTCAR.s
POTCAR.sb
```

