import scipy.constants
import ase.units


def kelvin_to_rydberg(temperature_K):
    k_B = scipy.constants.Boltzmann
    Ry_in_Joule = scipy.constants.Rydberg*scipy.constants.h*scipy.constants.c
    return (k_B*temperature_K)/Ry_in_Joule


def kelvin_to_eV(temperature_K):
    return kelvin_to_rydberg(temperature_K) * ase.units.Rydberg


def second_to_rydberg_time(time_s):
    Ry_in_Joule = scipy.constants.Rydberg*scipy.constants.h*scipy.constants.c
    t_0 = scipy.constants.hbar/Ry_in_Joule
    return time_s / t_0

