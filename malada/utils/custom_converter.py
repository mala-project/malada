"""Collection of custom converters."""

import scipy.constants
import ase.units


def kelvin_to_rydberg(temperature_K):
    """
    Convert a temperature from Kelvin to Rydberg energy units.

    Parameters
    ----------
    temperature_K : float
        Temperature in Kelvin.

    Returns
    -------
    temperature_Ry : float
        Temperature expressed in Rydberg.

    """
    k_B = scipy.constants.Boltzmann
    Ry_in_Joule = (
        scipy.constants.Rydberg * scipy.constants.h * scipy.constants.c
    )
    return (k_B * temperature_K) / Ry_in_Joule


def kelvin_to_eV(temperature_K):
    """
    Convert a temperature from Kelvin to electron volts.

    Parameters
    ----------
    temperature_K : float
        Temperature in Kelvin.

    Returns
    -------
    temperature_Ry : float
        Temperature expressed in electron volts.

    """
    return kelvin_to_rydberg(temperature_K) * ase.units.Rydberg


def second_to_rydberg_time(time_s):
    """
    Convert a time from second to Rydberg unit time.

    This is used in e.g. QE.

    Parameters
    ----------
    time_s : float
        Time in seconds.

    Returns
    -------
    time_ry : float
        Time in Rydberg unit time.

    """
    Ry_in_Joule = (
        scipy.constants.Rydberg * scipy.constants.h * scipy.constants.c
    )
    t_0 = scipy.constants.hbar / Ry_in_Joule
    return time_s / t_0
