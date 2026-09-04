"""Chemical model and agreement metrics."""

import numpy as np

MIN_DISTANCE = 200.0  # 2.00 Å, represented in cÅ coordinate units.


def delta_value(eu_hydrogen, eu_oxygen, oxygen_hydrogen):
    """Calculate the unscaled lanthanide-induced shift contribution."""
    cos_theta = (
        eu_hydrogen**2 + eu_oxygen**2 - oxygen_hydrogen**2
    ) / (2 * eu_hydrogen * eu_oxygen)
    return (3 * cos_theta**2 - 1) / eu_hydrogen**3


def r_factor(calc, experimental):
    """Calculate the agreement factor used by the reference model."""
    calc = np.asarray(calc, dtype=float)
    experimental = np.asarray(experimental, dtype=float)
    return 100.0 * np.sqrt(
        np.mean(((calc - experimental) / experimental) ** 2)
    )


def pearson_coefficient(calc, experimental):
    """Calculate Pearson correlation for one or many calculated vectors."""
    calc = np.asarray(calc, dtype=float)
    experimental = np.asarray(experimental, dtype=float)

    centered_calc = calc - calc.mean(axis=-1, keepdims=True)
    centered_exp = experimental - experimental.mean()

    numerator = np.sum(centered_calc * centered_exp, axis=-1)
    denominator = np.sqrt(
        np.sum(centered_calc**2, axis=-1)
        * np.sum(centered_exp**2)
    )

    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator != 0,
    )
