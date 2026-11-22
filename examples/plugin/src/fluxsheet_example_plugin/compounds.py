"""A compound the database lacks: a petroleum cut as a pseudo-component.

The definition has the same fields as a flowsheet's custom_compounds; what it does not give
(critical constants, heat capacity, vapor pressure) is estimated from the Joback groups.
"""

COMPOUNDS = [
    {
        "name": "light-naphtha-cut",
        "formula": "C8H18",
        "joback_groups": {"-CH3": 2, "-CH2-": 6},
        "constants": {"Tb": "120 degC"},
    },
]
