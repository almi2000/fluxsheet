"""Numerical comparison helpers; callers must separately qualify experimental data.

These functions do not run Fluxsheet or establish parameter independence.
"""

import math
from itertools import pairwise
from statistics import fmean, pstdev


def _values(values, name):
    result = [float(v) for v in values]
    if not result or not all(math.isfinite(v) for v in result):
        raise ValueError(f"{name} must contain finite values")
    return result


def interpolate(times, values, targets):
    """Linearly interpolate inside the supplied domain; never extrapolate."""
    times, values, targets = (_values(times, "times"), _values(values, "values"),
                              _values(targets, "targets"))
    if len(times) != len(values) or any(b <= a for a, b in pairwise(times)):
        raise ValueError("times must be strictly increasing and match values")
    if min(targets) < times[0] or max(targets) > times[-1]:
        raise ValueError("interpolation would extrapolate")
    from bisect import bisect_right
    out = []
    for t in targets:
        i = min(bisect_right(times, t) - 1, len(times) - 2)
        if len(times) == 1:
            out.append(values[0])
        else:
            w = (t - times[i]) / (times[i + 1] - times[i])
            out.append(values[i] + w * (values[i + 1] - values[i]))
    return out


def errors(measured, predicted, mean_limit, max_limit, relative=False, zero_floor=0.0):
    """Return pointwise and aggregate errors with fixed screening thresholds."""
    measured, predicted = _values(measured, "measured"), _values(predicted, "predicted")
    if not all(math.isfinite(x) and x >= 0 for x in (mean_limit, max_limit, zero_floor)):
        raise ValueError("limits must be finite and nonnegative")
    if len(measured) != len(predicted):
        raise ValueError("measurement and prediction lengths differ")
    if relative and any(abs(x) <= zero_floor for x in measured):
        raise ValueError("relative error is not checkable near zero")
    residuals = [p - m for m, p in zip(measured, predicted)]
    absolute = [abs(r / m) if relative else abs(r) for r, m in zip(residuals, measured)]
    maximum = max(absolute)
    mean = fmean(absolute)
    return {"measured": measured, "predicted": predicted, "residuals": residuals,
            "absolute_errors": absolute, "mean_absolute_error": mean,
            "maximum_absolute_error": maximum, "maximum_error_index": absolute.index(maximum),
            "passed": mean <= mean_limit and maximum <= max_limit}


def kinetics(measured, predicted):
    """Compare one run's dimensionless conversions using the registered limits."""
    return errors(measured, predicted, 0.03, 0.08)


def column(measured, predicted, quantity, *, zero_floor=0.0):
    """Compare one column quantity across operating points; duty errors are fractional."""
    limits = {"composition": (0.02, 0.05), "temperature": (2.0, 5.0),
              "duty": (0.1, 0.2)}
    if quantity not in limits:
        raise ValueError("unsupported column quantity")
    return errors(measured, predicted, *limits[quantity], relative=quantity == "duty",
                  zero_floor=zero_floor)


def dynamic_step(baseline_times, baseline, measured, predicted, end_start):
    """Evaluate one aligned step; end_start indexes the prespecified end window."""
    times = _values(baseline_times, "baseline times")
    base = _values(baseline, "baseline")
    measured = _values(measured, "measured")
    predicted = _values(predicted, "predicted")
    if len(base) < 3 or len(base) != len(times) or any(b <= a for a, b in pairwise(times)):
        raise ValueError("baseline needs at least three strictly ordered samples")
    if len(measured) != len(predicted) or not 0 <= end_start < len(measured):
        raise ValueError("invalid trajectory or end window")
    tmean, ymean = fmean(times), fmean(base)
    slope = sum((t-tmean)*(y-ymean) for t, y in zip(times, base)) / sum(
        (t-tmean)**2 for t in times)
    noise = pstdev([y - slope * (t-tmean) for t, y in zip(times, base)])
    amplitude = abs(fmean(measured[end_start:]) - ymean)
    if amplitude < 5 * noise or amplitude == 0:
        raise ValueError("step amplitude is less than five times baseline noise")
    residuals = [p-m for m, p in zip(measured, predicted)]
    rmse = math.sqrt(fmean([r*r for r in residuals]))
    end_error = abs(fmean(residuals[end_start:]))
    return {"measured": measured, "predicted": predicted, "residuals": residuals,
            "baseline_noise": noise, "step_amplitude": amplitude, "rmse": rmse,
            "end_error": end_error, "passed": max(rmse, end_error) <= 0.1 * amplitude}
