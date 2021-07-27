import scipy.constants

def kelvin_to_rydberg(temperature_K):
    k_B = scipy.constants.Boltzmann
    Ry_in_Joule = scipy.constants.Rydberg*scipy.constants.h*scipy.constants.c
    return (k_B*temperature_K)/Ry_in_Joule
