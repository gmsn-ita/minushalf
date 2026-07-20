"""
Parse the `projwfc.x` output file from Quantum ESPRESSO.

Read and process the spin-up projection file (`projwfc_up`), which
contains the projections of Kohn-Sham states onto atomic orbitals as
computed by `projwfc.x`.
"""
import re
from collections import defaultdict
from minushalf.softwares.band_projection_file import BandProjectionFile


class ProjOutput(BandProjectionFile):
    """
    Parse a `projwfc_up` file and store the band projection information.
    File structure
    --------------
    The file consists of blocks, one per atomic state, each introduced by a
    header line with the following fields:
        state_idx  atom_idx  symbol  wfc_label  wfc_idx  l  m
    Example::
        1    1 Al   3S     1    0    1
        2    1 Al   3P     2    1    1
    Each header is followed by rows of k-point and band projections:
        kpoint_index    band_index    projection_value
    Example::
        1       1        0.4379301244
        1       2        0.1468167336
        ...
        40      20       0.0003262801
    The projection value is :math:`|\\langle \\psi_{nk} | \\phi_i \\rangle|^2`
    (already squared).
    """

    # ------------------------------------------------------------------ #
    #  Header line of a state block:                                       #
    #  "   1    1 Al   3S     1    0    1"                                 #
    #  groups: state_idx  atom_idx  symbol  wfc_label  wfc_idx  l  m      #
    # ------------------------------------------------------------------ #
    _STATE_HEADER_REGEX = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+([A-Za-z]+)\s+\S+\s+\d+\s+(\d+)\s+(\d+)\s*$"
    )

    # Data row inside a state block:  kpoint  band  value
    _DATA_ROW_REGEX = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)\s*$"
    )

    def __init__(self, filename: str):
        """
        Args:
            filename (str): path to the projwfc_up (or projwfc_down) file
        Members:
            filename     : stored path
            num_kpoints  : number of k-points found in the file
            num_bands    : number of bands found in the file
            num_states   : total number of atomic projection states
            states_info  : list of dicts, one per state, with keys
                               state_idx, atom_idx, symbol, l, m
            _projections : nested dict  [state_idx][kpoint][band] = value
                           populated lazily on first call to get_band_projection
        """
        self.filename = filename
        self.states_info = []
        self._projections = {}

        self.num_kpoints, self.num_bands, self.num_states = (
            self._get_dimensions()
        )

    # ------------------------------------------------------------------ #
    #  Public interface (mirrors Procar.get_band_projection)               #
    # ------------------------------------------------------------------ #

    def get_band_projection(self, kpoint: int, band_number: int) -> dict:
        """
        Return the projection of a given (kpoint, band) onto every
        atomic state.

        Mirrors the Procar interface:  the return value is a dict whose
        keys are atom indices (str) and whose values are lists of floats,
        one entry per orbital of that atom, ordered by increasing m.

        Args:
            kpoint      (int): 1-based k-point index
            band_number (int): 1-based band index

        Returns:
            projections (dict):
                { "1": [proj_state1, proj_state2, ...],   # atom 1 orbitals
                  "2": [proj_state5, proj_state6, ...],   # atom 2 orbitals
                  ... }
        """
        if not self._projections:
            self._load_projections()

        projections = {}
        for state in self.states_info:
            atom_key = str(state["atom_idx"])
            value = (
                self._projections
                .get(state["state_idx"], {})
                .get(kpoint, {})
                .get(band_number, 0.0)
            )
            projections.setdefault(atom_key, []).append(value)

        return projections

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _get_dimensions(self) -> tuple:
        """
        Single-pass scan to determine:
          - number of k-points  (from header line: "nstates  nkpts  nbands")
          - number of bands     (from header line: "nstates  nkpts  nbands")
          - number of states    (count of state header lines)

        The file header contains a summary line of exactly three integers:
            16      40      20
        meaning 16 atomic states, 40 k-points, 20 bands. This is read
        directly rather than inferring from max indices in the data rows,
        which would risk picking up stray integers from the file header.

        Also populates self.states_info as a side-effect so we only
        read the file once during __init__.

        Returns:
            (num_kpoints, num_bands, num_states) (tuple[int, int, int])
        """
        # Matches exactly three integers on a line — the summary header
        summary_regex = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s*$")

        num_kpoints = None
        num_bands   = None
        states_seen = set()

        with open(self.filename, "r") as fh:
            for line in fh:
                # Parse summary header before state blocks begin
                if num_kpoints is None:
                    summary_match = summary_regex.match(line)
                    if summary_match:
                        # format: num_states  num_kpoints  num_bands
                        num_kpoints = int(summary_match.group(2))
                        num_bands   = int(summary_match.group(3))
                    continue

                # Once header is found, collect state header lines
                state_match = self._STATE_HEADER_REGEX.match(line)
                if state_match:
                    state_idx = int(state_match.group(1))
                    if state_idx not in states_seen:
                        states_seen.add(state_idx)
                        self.states_info.append({
                            "state_idx": state_idx,
                            "atom_idx":  int(state_match.group(2)),
                            "symbol":    state_match.group(3),
                            "l":         int(state_match.group(4)),
                            "m":         int(state_match.group(5)),
                        })

        if num_kpoints is None or num_bands is None:
            raise Exception(
                "ProjOutput parser could not find the summary header line "
                f"(nstates nkpts nbands) in {self.filename}"
            )
        if not states_seen:
            raise Exception(
                "ProjOutput parser could not find any state headers "
                f"in {self.filename}"
            )

        return num_kpoints, num_bands, len(states_seen)

    def _load_projections(self) -> None:
        """
        Full single-pass load of all projection values into
        self._projections[state_idx][kpoint][band] = value.

        Called lazily on the first get_band_projection call so that
        __init__ stays cheap (only _get_dimensions runs at construction).
        """
        # nested defaultdict for ergonomic assignment
        raw = defaultdict(lambda: defaultdict(dict))

        current_state = None
        with open(self.filename, "r") as fh:
            for line in fh:
                state_match = self._STATE_HEADER_REGEX.match(line)
                if state_match:
                    current_state = int(state_match.group(1))
                    continue

                if current_state is None:
                    continue

                data_match = self._DATA_ROW_REGEX.match(line)
                if data_match:
                    kpt   = int(data_match.group(1))
                    band  = int(data_match.group(2))
                    value = float(data_match.group(3))
                    raw[current_state][kpt][band] = value

        # convert to plain dicts for a stable, serialisable structure
        self._projections = {
            s: {k: dict(bands) for k, bands in kpts.items()}
            for s, kpts in raw.items()
        }
        