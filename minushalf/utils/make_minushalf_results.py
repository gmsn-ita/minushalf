"""
function make minushalf_results.dat
"""
from loguru import logger


def make_minushalf_results(
    gap: float,
    valence_cuts: dict,
    conduction_cuts: dict = None,
    name: str = "minushalf_results.dat",
) -> None:
    """
    Make output file for execute command, minushalf_results.dat.

        Args:
            gap (float): final gap in the correction method

            valence_cuts (dict): dictionary inform the atom symbol, orbital
            and cut for valence correction in the following format.
            {
                (symbol,orbital):cut
            }

            conduction_cuts (dict):dictionary inform the atom symbol, orbital
            and cut for condsuction correction in the following format.
            {
                (symbol,orbital):cut
            }

            name (str): name of the file
    """
    logger.info("Writing output file")
    with open(name, "w") as file:
        ## Write valence cuts
        file.write("Valence correction cuts:\n")
        for key, value in valence_cuts.items():
            symbol, orbital = key
            cut = value
            file.write("\t({},{}):{:.2f}A\n".format(symbol, orbital, cut))
        file.write(
            "----------------------------------------------------------------\n"
        )
        ## Write conduction cuts
        if conduction_cuts:
            file.write("Conduction correction cuts:\n")
            for key, value in conduction_cuts.items():
                symbol, orbital = key
                cut = value
                file.write("\t({},{}):{:.2f}A\n".format(symbol, orbital, cut))
            file.write(
                "----------------------------------------------------------------\n"
            )
        ## Write gap
        file.write("GAP: {:.3}eV".format(gap))
