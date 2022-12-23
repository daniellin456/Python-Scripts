import numpy as np


def CII_Cooling(nH, T, x):
    CII = 92.0 * 1.38e-16 * 2.0 * (
            (2.8E-7 * (T / 100.0) ** -0.5) * x + 8e-10 * (T / 100.0) ** 0.07) * 3.5e-4 * 0.4 * np.exp(-92.0 / T)
    return CII


def OI_Cooling(nH, T, x):
    OI = 1e-26 * np.sqrt(T) * (24.0 * np.exp(-228.0 / T) + 7.0 * np.exp(-326.0 / T)) * 4.5E-4
    return OI


def H2_Cooling(nH, T, x):
    H2 = 7.3E-19 * x * np.exp(-118400.0 / T)
    return H2


def OI_Line_Cooling(nH, T, x):
    OI_Line = 4.5e-4 * (2.3e4 * 1.38e-16 / 3.0 * (5.1e-9 * (T / 10000.) ** 0.57 * x + 1e-12) * np.exp(
        -2.3e4 / T) + 4.9e4 * 1.38e-16 / 3.0 * (2.5e-9 * (T / 10000.) ** 0.57 * x + 1e-12) * np.exp(
        -4.9e4 / T) + 2.6e4 * 1.38e-16 * 1.0 * (5.2e-9 * (T / 10000.) ** 0.57 * x + 1e-12) * np.exp(
        -2.6e4 / T))

    return OI_Line


def CII_Line_Cooling(nH, T, x):
    CII_Line = 6.2e4 * 1.38e-16 * 1.0 * (2.3e-8 * (T / 10000.0) ** -0.5 * x + 1e-12) * np.exp(
        -6.2e4 / T) * 3.5e-4 * 0.4
    return CII_Line


def Recombination_Cooling(nH, T, x, param, bet):
    Recombination = 4.65E-30 * (T ** 0.94) * (param ** bet) * x
    return Recombination


def UV_Heating(nH, T, x, param, G0):
    epsilon = 4.9E-2 / (1.0 + (param / 1925.) ** 0.73) + 3.7E-2 * (T / 1e4) ** 0.7 / (1. + (param / 5e3))
    return 1e-24 * epsilon * G0


def cooling(nH, T, x):
    G0 = 1.0 / 1.7
    bet = 0.74 / (T ** 0.068)
    param = G0 * np.sqrt(T) / (nH * x)
    Total_Cooling = CII_Cooling(nH, T, x) + OI_Cooling(nH, T, x) + H2_Cooling(nH, T, x) + OI_Line_Cooling(nH, T, x) + \
                    CII_Line_Cooling(nH, T, x) + Recombination_Cooling(nH, T, x, param, bet)
    return Total_Cooling


def heating(nH, T, x):
    G0 = 1.0 / 1.7
    param = G0 * np.sqrt(T) / (nH * x)
    return UV_Heating(nH, T, x, param, G0)


def dCoolingdRho(nH, T, x):
    bet = 0.74 / (T ** 0.068)
    if x >= 0.1:
        dCII_CoolingdRho = 0
        dOI_CoolingdRho = 0
        dH2_CoolingdRho = 0
        dCII_Line_CoolingdRho = 0
        dOI_Line_CoolingdRho = 0
        dRecombinationdRho = -3.441e-31 * T ** 0.872 * (5.88235294117647 * T ** 0.5 / nH) ** bet / nH
    elif x <= 1.4e-4:
        dCII_CoolingdRho = 0
        dOI_CoolingdRho = 0
        dH2_CoolingdRho = 0
        dCII_Line_CoolingdRho = 0
        dOI_Line_CoolingdRho = 0
        dRecombinationdRho = -4.8174e-34 * T ** 0.872 * (4201.68067226891 * T ** 0.5 / nH) ** bet / nH
    else:
        dCII_CoolingdRho = -1.51085996659307e-26 * np.exp(-92.0 / T) / (T ** 0.25 * nH ** 2)
        dOI_CoolingdRho = 0
        dH2_CoolingdRho = -1.108062092123e-21 * T ** 0.25 * np.exp(-118400.0 / T) / nH ** 2
        dCII_Line_CoolingdRho = -4.18184455039153e-24 * np.exp(-6.2e4) / (T ** 0.25 * nH ** 2)
        dOI_Line_CoolingdRho = -1.93423659159274e-29 * T ** 0.82 * np.exp(
            -2.3e4 / T) * nH ** -2 - 2.0199828002567e-29 * T ** 0.82 * np.exp(
            -4.9e4 / T) * nH ** -2 - 6.68820427578872e-29 * T ** 0.82 * np.exp(-2.6e4 / T) * nH ** -2
        dRecombinationdRho = -7.05820373749582e-33 * T ** 1.19 * (387.534026981419 * T ** 0.25) ** bet / nH ** 2

    dTotal_CoolingdRho = dCII_CoolingdRho + dOI_CoolingdRho + dH2_CoolingdRho + dCII_Line_CoolingdRho + dOI_Line_CoolingdRho + dRecombinationdRho
    return dTotal_CoolingdRho


def dHeatingdRho(nH, T, x):
    if x >= 0.1:
        dUV_HeatingdRho = 4.05820402221877e-32 * T ** 1.2 / (
                nH ** 2 * (0.00117647058823529 * T ** 0.5 / nH + 1.0) ** 2) + 3.07049129489545e-28 * (
                                  T ** 0.5 / nH) ** 0.73 / (
                                  nH * (0.0145927738365174 * (T ** 0.5 / nH) ** 0.73 + 1.0) ** 2)
    elif x <= 1.4e-4:
        dUV_HeatingdRho = 2.89871715872769e-29 * T ** 1.2 / (
                nH ** 2 * (0.840336134453782 * T ** 0.5 / nH + 1.0) ** 2) + 1.19016077589066e-26 * (
                                  T ** 0.5 / nH) ** 0.73 / (nH * ((T ** 0.5 / nH) ** 0.73 + 0.565634140065452) ** 2)
    else:
        dUV_HeatingdRho = 0
    return dUV_HeatingdRho


def dCoolingdTemp(nH, T, x):
    bet = 0.74 / (T ** 0.068)
    if x >= 0.1:
        dCII_CoolingdTemp = (-4.976832e-25 / T ** 1.5 + 1.44215842722077e-28 / T ** 0.93) * np.exp(-92.0 / T) + 92.0 * (
                9.953664e-25 / T ** 0.5 + 2.0602263246011e-27 * T ** 0.07) * np.exp(-92.0 / T) / T ** 2
        dOI_CoolingdTemp = 2.25e-30 * (
                7.0 * np.exp(-326.0 / T) + 24.0 * np.exp(-228.0 / T)) / T ** 0.5 + 4.5e-30 * T ** 0.5 * (
                                   2282.0 * np.exp(-326.0 / T) / T ** 2 + 5472.0 * np.exp(-228.0 / T) / T ** 2)
        dH2_CoolingdTemp = 8.6432e-15 * np.exp(-118400.0 / T) / T ** 2
        dCII_Line_CoolingdTemp = -1.377516e-22 * np.exp(-62000.0 / T) / T ** 1.5 + 8.68 * (
                1.96788e-18 / T ** 0.5 + 8.556e-24) * np.exp(-62000.0 / T) / T ** 2
        dOI_Line_CoolingdTemp = 7.58544894877166e-28 * np.exp(-49000.0 / T) / T ** 0.43 + 2.51155762744636e-27 * np.exp(
            -26000.0 / T) / T ** 0.43 + 7.26345438115033e-28 * np.exp(-23000.0 / T) / T ** 0.43 + 10.35 * (
                                        2.83175609401572e-24 * T ** 0.57 + 1.058e-24) * np.exp(
            -23000.0 / T) / T ** 2 + 22.05 * (2.95729003850747e-24 * T ** 0.57 + 2.254e-24) * np.exp(
            -49000.0 / T) / T ** 2 + 11.7 * (9.79164767035616e-24 * T ** 0.57 + 3.588e-24) * np.exp(
            -26000.0 / T) / T ** 2
        dRecombination_CoolingdTemp = 4.371e-31 * (
                5.88235294117647 * T ** 0.5 / nH) ** bet / T ** 0.0600000000000001 + 4.65e-31 * T ** 0.94 * (
                                              5.88235294117647 * T ** 0.5 / nH) ** bet * (-0.05032 * np.log(
            5.88235294117647 * T ** 0.5 / nH) / T ** 1.068 + 0.37 / T ** 1.068)
    elif x <= 1.4e-4:
        dCII_CoolingdTemp = (-6.9675648e-28 / T ** 1.5 + 1.44215842722077e-28 / T ** 0.93) * np.exp(
            -92.0 / T) + 92.0 * (1.39351296e-27 / T ** 0.5 + 2.0602263246011e-27 * T ** 0.07) * np.exp(
            -92.0 / T) / T ** 2
        dOI_CoolingdTemp = 2.25e-30 * (
                7.0 * np.exp(-326.0 / T) + 24.0 * np.exp(-228.0 / T)) / T ** 0.5 + 4.5e-30 * T ** 0.5 * (
                                   2282.0 * np.exp(-326.0 / T) / T ** 2 + 5472.0 * np.exp(-228.0 / T) / T ** 2)
        dH2_CoolingdTemp = 1.210048e-17 * np.exp(-118400.0 / T) / T ** 2
        dCII_Line_CoolingdTemp = -1.9285224e-25 * np.exp(-62000.0 / T) / T ** 1.5 + 8.68 * (
                2.755032e-21 / T ** 0.5 + 8.556e-24) * np.exp(-62000.0 / T) / T ** 2
        dOI_Line_CoolingdTemp = 1.06196285282803e-30 * np.exp(-49000.0 / T) / T ** 0.43 + 3.5161806784249e-30 * np.exp(
            -26000.0 / T) / T ** 0.43 + 1.01688361336105e-30 * np.exp(-23000.0 / T) / T ** 0.43 + 10.35 * (
                                        3.96445853162201e-27 * T ** 0.57 + 1.058e-24) * np.exp(
            -23000.0 / T) / T ** 2 + 22.05 * (4.14020605391046e-27 * T ** 0.57 + 2.254e-24) * np.exp(
            -49000.0 / T) / T ** 2 + 11.7 * (1.37083067384986e-26 * T ** 0.57 + 3.588e-24) * np.exp(
            -26000.0 / T) / T ** 2
        dRecombination_CoolingdTemp = 6.1194e-34 * (
                4201.68067226891 * T ** 0.5 / nH) ** bet / T ** 0.0600000000000001 + 6.51e-34 * T ** 0.94 * (
                                              4201.68067226891 * T ** 0.5 / nH) ** bet * (-0.05032 * np.log(
            4201.68067226891 * T ** 0.5 / nH) / T ** 1.068 + 0.37 / T ** 1.068)
    else:
        dCII_CoolingdTemp = (-3.77714991648267e-27 / (T ** 1.25 * nH) + 1.44215842722077e-28 / T ** 0.93) * np.exp(
            -92.0 / T) + 92.0 * (1.51085996659307e-26 / (T ** 0.25 * nH) + 2.0602263246011e-27 * T ** 0.07) * np.exp(
            -92.0 / T) / T ** 2
        dOI_CoolingdTemp = 2.25e-30 * (
                7.0 * np.exp(-326.0 / T) + 24.0 * np.exp(-228.0 / T)) / T ** 0.5 + 4.5e-30 * T ** 0.5 * (
                                   2282.0 * np.exp(-326.0 / T) / T ** 2 + 5472.0 * np.exp(-228.0 / T) / T ** 2)
        dH2_CoolingdTemp = 1.31194551707363e-16 * np.exp(-118400.0 / T) / (
                T ** 1.75 * nH) + 2.7701552303075e-22 * np.exp(-118400.0 / T) / (T ** 0.75 * nH)
        dCII_Line_CoolingdTemp = -1.04546113759788e-24 * np.exp(-62000.0 / T) / (T ** 1.25 * nH) + 8.68 * (
                2.98703182170823e-20 / (T ** 0.25 * nH) + 8.556e-24) * np.exp(-62000.0 / T) / T ** 2
        dOI_Line_CoolingdTemp = 1.6563858962105e-29 * np.exp(-49000.0 / T) / (
                T ** 0.18 * nH) + 5.48432750614675e-29 * np.exp(-26000.0 / T) / (
                                        T ** 0.18 * nH) + 1.58607400510605e-29 * np.exp(-23000.0 / T) / (
                                        T ** 0.18 * nH) + 10.35 * (
                                        4.29830353687276e-26 * T ** 0.82 / nH + 1.058e-24) * np.exp(
            -23000.0 / T) / T ** 2 + 22.05 * (4.48885066723712e-26 * T ** 0.82 / nH + 2.254e-24) * np.exp(
            -49000.0 / T) / T ** 2 + 11.7 * (1.48626761684194e-25 * T ** 0.82 / nH + 3.588e-24) * np.exp(
            -26000.0 / T) / T ** 2
        dRecombination_CoolingdTemp = 8.39926244762003e-33 * T ** 0.19 * (
                387.534026981419 * T ** 0.25) ** bet / nH + 7.05820373749582e-33 * T ** 1.19 * (
                                              387.534026981419 * T ** 0.25) ** bet * (-0.05032 * np.log(
            387.534026981419 * T ** 0.25) / T ** 1.068 + 0.185 / T ** 1.068) / nH

    dTotal_CoolingdTemp = dCII_CoolingdTemp + dOI_CoolingdTemp + dH2_CoolingdTemp + dCII_Line_CoolingdTemp + dOI_Line_CoolingdTemp + dRecombination_CoolingdTemp
    return dTotal_CoolingdTemp


def dHeatingdTemp(nH, T, x):
    if x >= 0.1:
        dUV_HeatingdTemp = -1.53524564744772e-28 * (T ** 0.5 / nH) ** 0.73 / (
                T ** 1.0 * (0.0145927738365174 * (T ** 0.5 / nH) ** 0.73 + 1.0) ** 2) + 2.41463139322017e-29 / (
                                   T ** 0.3 * (
                                   0.00117647058823529 * T ** 0.5 / nH + 1.0)) - 2.02910201110938e-32 * T ** 0.2 / (
                                   nH * (0.00117647058823529 * T ** 0.5 / nH + 1.0) ** 2)
    elif x <= 1.4e-4:
        dUV_HeatingdTemp = -5.9508038794533e-27 * (T ** 0.5 / nH) ** 0.73 / (
                T ** 1.0 * ((T ** 0.5 / nH) ** 0.73 + 0.565634140065452) ** 2) + 2.41463139322017e-29 / (
                                   T ** 0.3 * (
                                   0.840336134453782 * T ** 0.5 / nH + 1.0)) - 1.44935857936385e-29 * T ** 0.2 / (
                                   nH * (0.840336134453782 * T ** 0.5 / nH + 1.0) ** 2)
    else:
        dUV_HeatingdTemp = -1.63245709529796e-27 / (
                T ** 0.8175 * (0.310335707241436 * T ** 0.1825 + 1.0) ** 2) + 2.41463139322017e-29 / (
                                   T ** 0.3 * (0.0775068053962838 * T ** 0.25 + 1.0)) - 6.68394162493118e-31 / (
                                   T ** 0.05 * (0.0775068053962838 * T ** 0.25 + 1.0) ** 2)
    return dUV_HeatingdTemp
