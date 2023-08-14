import pathlib

import cobra
import cobra.manipulation



def ecoli_dh5a(variant: str = "base") -> cobra.Model:

    #root = pathlib.Path(__file__).parents[2] / "data"

    if variant in {"base", "corrected"}:
        sbml_path = "iEC1368_DH5a.xml"
        return cobra.io.read_sbml_model(str(sbml_path))

    else:
        raise ValueError()
