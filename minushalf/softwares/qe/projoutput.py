"""
Reads projwfc_up file, an output of
Quantum ESPRESSO projwfc.x software
"""
import re
from itertools import islice
from collections import defaultdict
from minushalf.softwares.band_projection_file import BandProjectionFile


class ProjOutput(BandProjectionFile):
    """
    Reads a projwfc_up file and stores band projection information.

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
    #  Private helpers                                                   #
    # ------------------------------------------------------------------ #

    def _get_dimensions(self) -> tuple:
        """
        Single-pass scan to determine:
          - number of k-points  (max kpoint index seen)
          - number of bands     (max band index seen)
          - number of states    (count of state header lines)

        Also populates self.states_info as a side-effect so we only
        read the file once during __init__.

        Returns:
            (num_kpoints, num_bands, num_states) (tuple[int, int, int])
        """
        max_kpoint = 0
        max_band = 0
        states_seen = set()

        with open(self.filename, "r") as fh:
            for line in fh:
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
                    continue

                data_match = self._DATA_ROW_REGEX.match(line)
                if data_match:
                    kpt  = int(data_match.group(1))
                    band = int(data_match.group(2))
                    if kpt  > max_kpoint:
                        max_kpoint = kpt
                    if band > max_band:
                        max_band = band

        if not states_seen:
            raise Exception(
                "ProjOutput parser could not find any state headers "
                f"in {self.filename}"
            )
        if max_kpoint == 0 or max_band == 0:
            raise Exception(
                "ProjOutput parser could not find any projection data "
                f"in {self.filename}"
            )

        return max_kpoint, max_band, len(states_seen)

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