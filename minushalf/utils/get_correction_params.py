"""
Extract the parameters for the correction
"""
import pandas as pd
from loguru import logger
from minushalf.io.minushalf_yaml_interface import MinushalfYaml
from minushalf.softwares.software_abstract_factory import SoftwaresAbstractFactory
from minushalf.utils.cut_initial_guess import CutInitialGuess
from minushalf.utils.fractionary_correction_indexes import get_fractionary_correction_indexes
from minushalf.utils.simple_correction_indexes import get_simple_correction_indexes
from minushalf.utils.cut_initial_guess_methods import (CutInitialGuessMethods)
from minushalf.utils.band_structure import BandStructure
from minushalf.utils.projection_to_df import projection_to_df


def _get_vbm_projection(factory: SoftwaresAbstractFactory, is_indirect: bool) -> pd.DataFrame:
    """
    Returns vbm projection
    """
    band_structure = BandStructure.create(factory)
    vbm_projection = band_structure.vbm_projection(is_indirect=is_indirect)
    normalized_df = projection_to_df(vbm_projection)
    return normalized_df


def _get_cbm_projection(factory: SoftwaresAbstractFactory, is_indirect: bool) -> pd.DataFrame:
    """
    Returns cbm projection
    """
    band_structure = BandStructure.create(factory)
    cbm_projection = band_structure.cbm_projection(is_indirect=is_indirect)
    normalized_df = projection_to_df(cbm_projection)
    return normalized_df


def _overwrite_band_projection(new_values: list,
                               band_projection: pd.DataFrame) -> pd.DataFrame:
    """
    Overwrite values in VBM or CBM band projection
        Args:
            new_values (list): Arguments passed in overwrite_vbm or overwrite_cbm.
            band_projection (pd.Dataframe): Dataframe with the value of projections
                                            of VBM or CBM.
        Returns:
            modified_band_projection (pd.Dataframe): Band projection data frame
                                                    with the values overwrited.
    """
    for case in new_values:
        atom = case[0].capitalize()
        orbital = case[1].lower()
        projection = int(case[2])
        band_projection[orbital][atom] = projection
    return band_projection


def _get_band_characters(band_location: list,
                         factory: SoftwaresAbstractFactory) -> pd.DataFrame:
    """
    Overwrite values in VBM or CBM band projection

        Args:

            band_location (List[int]): Kpoint and band number, respectively.

            factory (SoftwaresAbstractFactory): Software factory to create band structure class

        Returns:

            modified_band_projection (pd.Dataframe): Band projection data frame
                                                    with the values overwrited.

    """

    band_structure = BandStructure.create(factory)
    vbm_projection = band_structure.band_projection(*band_location)
    normalized_df = projection_to_df(vbm_projection)
    return normalized_df


def _get_valence_band_projection(minushalf_yaml: MinushalfYaml,
                                 software_factory: SoftwaresAbstractFactory):
    """
    Implements logic to return valence band projection
    """
    # Specify that the VBM should be another band
    has_overwrite = bool(minushalf_yaml.get_overwrite_vbm())

    projection_df = None
    if has_overwrite:
        band_location = minushalf_yaml.get_overwrite_vbm()
        projection_df = _get_band_characters(band_location, software_factory)
    else:
        projection_df = _get_vbm_projection(
            software_factory, is_indirect=minushalf_yaml.get_indirect())
    # If the vbm characters are overwritten manually
    has_replace = bool(minushalf_yaml.get_vbm_characters())
    if has_replace:
        projection_df = _overwrite_band_projection(
            minushalf_yaml.get_vbm_characters(), projection_df)

    return projection_df


def _get_conduction_band_projection(
        minushalf_yaml: MinushalfYaml,
        software_factory: SoftwaresAbstractFactory):
    """
    Implements logic to return conduction band projection
    """
    has_overwrite = bool(minushalf_yaml.get_overwrite_cbm())

    projection_df = None
    if has_overwrite:
        band_location = minushalf_yaml.get_overwrite_cbm()
        projection_df = _get_band_characters(band_location, software_factory)
    else:
        projection_df = _get_cbm_projection(
            software_factory, is_indirect=minushalf_yaml.get_indirect())

    has_replace = bool(minushalf_yaml.get_cbm_characters())
    if has_replace:
        projection_df = _overwrite_band_projection(
            minushalf_yaml.get_cbm_characters(), projection_df)

    return projection_df


def _guess_distance(symbol, software_factory):
    """
    Guess the distance
    """
    atoms_map = software_factory.get_atoms_map()
    ion_index = None
    for key, value in atoms_map.items():
        if value == symbol:
            ion_index = key
            break

    nearest_distance = software_factory.get_nearest_neighbor_distance(
        ion_index)
    cut_guesser = CutInitialGuess()
    return cut_guesser.guess(nearest_distance,
                             CutInitialGuessMethods.three_dimensions.value)


def _get_cut_initial_guess(initial_guess: list, correction_indexes: dict,
                           software_factory: SoftwaresAbstractFactory) -> dict:
    """
    Get cut initial guess parameter
    """
    if not initial_guess:
        initial_guess = []
    cut_guesses = {(e[0], e[1]): e[2] for e in initial_guess}
    for atom, orbitals in correction_indexes.items():
        for orbital in orbitals:
            if not (atom, orbital) in cut_guesses:
                cut_guesses[(atom, orbital)] = _guess_distance(
                    atom, software_factory)
    return cut_guesses


def _get_divide_character(divide_characters: list, correction_indexes: dict,
                          software_factory: SoftwaresAbstractFactory):
    """
    Get divide character param
    """
    if not divide_characters:
        divide_characters = []

    dividers = {(e[0], e[1]): e[2] for e in divide_characters}
    atoms_map = software_factory.get_atoms_map()
    for atom, orbitals in correction_indexes.items():
        for orbital in orbitals:
            if not (atom, orbital) in dividers:
                dividers[(
                    atom,
                    orbital)] = software_factory.get_number_of_equal_neighbors(
                        atoms_map=atoms_map, symbol=atom)
    return dividers


def _get_valence_correction_indexes(correction_code, band_projection, treshold: int):
    """
    Get the correction indexes needed to the corrrection
    """
    if "vf" in correction_code:
        return get_fractionary_correction_indexes(band_projection, treshold)
    else:
        return get_simple_correction_indexes(band_projection)


def _get_conduction_correction_indexes(correction_code, band_projection, treshold: int):
    """
    Get the correction indexes needed to the corrrection
    """
    if "cf" in correction_code:
        return get_fractionary_correction_indexes(band_projection, treshold)
    else:
        return get_simple_correction_indexes(band_projection)


def get_valence_correction_params(
    minushalf_yaml: MinushalfYaml,
    software_factory: SoftwaresAbstractFactory,
    **kwargs,
):
    """
    Returns the parameters for the valence correction
    """
    correction_code = minushalf_yaml.get_correction_code()
    params = kwargs
    params["software_factory"] = software_factory
    params["potential_filename"] = software_factory.get_potential_class(
    ).get_name()
    params["band_projection"] = _get_valence_band_projection(
        minushalf_yaml, software_factory)
    params["potential_folder"] = minushalf_yaml.get_potential_folder()
    params[
        "exchange_correlation_type"] = minushalf_yaml.get_exchange_corr_code()
    params["max_iterations"] = minushalf_yaml.get_max_iterations()
    params["calculation_code"] = minushalf_yaml.get_calculation_code()
    params["amplitude"] = minushalf_yaml.get_amplitude()
    params["tolerance"] = minushalf_yaml.get_tolerance()
    params["corrected_potfiles_folder"] = ".minushalf/corrected_valence_potfiles"
    params["correction_type"] = "valence"
    params["is_conduction"] = False
    params["indirect"] = minushalf_yaml.get_indirect()
    params["input_files"] = [
        "INCAR", "POSCAR", "KPOINTS", "CHGCAR"
    ] if params["indirect"] else ["INCAR", "POSCAR", "KPOINTS"]

    params["correction_indexes"] = _get_valence_correction_indexes(
        correction_code, params["band_projection"], minushalf_yaml.get_fractional_valence_treshold())
    params["cut_initial_guess"] = _get_cut_initial_guess(
        minushalf_yaml.get_valence_cut_initial_guess(),
        params["correction_indexes"], software_factory)

    params["divide_character"] = _get_divide_character(
        minushalf_yaml.get_divide_character(), params["correction_indexes"],
        software_factory)

    return params


def get_conduction_correction_params(
    minushalf_yaml: MinushalfYaml,
    software_factory: SoftwaresAbstractFactory,
    **kwargs,
):
    """
    Returns the parameters for the conduction correction
    """
    correction_code = minushalf_yaml.get_correction_code()
    params = kwargs
    params["software_factory"] = software_factory
    params["potential_filename"] = software_factory.get_potential_class(
    ).get_name()
    params["band_projection"] = _get_conduction_band_projection(
        minushalf_yaml, software_factory)
    params["potential_folder"] = ".minushalf/corrected_valence_potfiles"
    params[
        "exchange_correlation_type"] = minushalf_yaml.get_exchange_corr_code()
    params["max_iterations"] = minushalf_yaml.get_max_iterations()
    params["calculation_code"] = minushalf_yaml.get_calculation_code()
    params["amplitude"] = minushalf_yaml.get_amplitude()
    params["tolerance"] = minushalf_yaml.get_tolerance()
    params["corrected_potfiles_folder"] = ".minushalf/corrected_conduction_potfiles"
    params["correction_type"] = "conduction"
    params["is_conduction"] = True
    params["indirect"] = minushalf_yaml.get_indirect()
    params["input_files"] = [
        "INCAR", "POSCAR", "KPOINTS", "CHGCAR"
    ] if params["indirect"] else ["INCAR", "POSCAR", "KPOINTS"]
    params["correction_indexes"] = _get_conduction_correction_indexes(
        correction_code, params["band_projection"], minushalf_yaml.get_fractional_conduction_treshold())
    params["cut_initial_guess"] = _get_cut_initial_guess(
        minushalf_yaml.get_conduction_cut_initial_guess(),
        params["correction_indexes"], software_factory)

    params["divide_character"] = _get_divide_character(
        minushalf_yaml.get_divide_character(), params["correction_indexes"],
        software_factory)

    return params
