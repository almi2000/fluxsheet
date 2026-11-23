"""A property method: Hildebrand's regular-solution model for non-polar liquids."""

from thermo.regular_solution import RegularSolution

from fluxsheet.plugins import ActivityMethod


def regular_solution(consts, components):
    """Activity coefficients from the solubility parameters and liquid molar volumes of the
    compound database (thermo's ChemicalConstantsPackage), for the compounds at the indices
    ``components``: the solvents, since dissolved light gases (hydrogen, methane, ...) are
    Henry components (``henry=True`` below) and need neither."""
    missing = [consts.names[i] for i in components
               if not consts.solubility_parameters[i] or not consts.Vml_STPs[i]]
    if missing:
        raise ValueError(f"no solubility parameter or liquid volume for {', '.join(missing)}")
    n = len(components)
    return RegularSolution(xs=[1.0 / n] * n,
                           Vs=[consts.Vml_STPs[i] for i in components],
                           SPs=[consts.solubility_parameters[i] for i in components])


REGULAR_SOLUTION = ActivityMethod(
    key="regular-solution",
    label="Regular solution (Hildebrand)",
    description="Activity coefficients from solubility parameters and liquid volumes, with "
                "an ideal vapor. For non-polar mixtures at low pressure; no binary parameters. "
                "Light gases dissolve by Henry's law. From the example plug-in.",
    gibbs_excess=regular_solution,
    henry=True,
)
