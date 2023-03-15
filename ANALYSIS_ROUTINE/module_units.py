from pymses.utils import constants as C


def set_units(data):
    units = dict()
    print('====== Output Info ======')
    for key in data.info:
        if key == "dom_decomp_Hilbert_keys":
            continue
        print(key, data.info[key])
    print
    ""
    boxlen = data.info["boxlen"]
    mu = 2.31  # 2.31  #ramses2.31, pymses 1/0.76
    fac_mu = 2.31 * 0.76
    fac_len = boxlen / data.info["unit_length"].express(C.pc)
    data.info["unit_length"] = data.info["unit_length"] * fac_len
    data.info["unit_density"] = data.info["unit_density"]
    data.info["unit_mass"] = data.info['unit_density'] * \
                             data.info['unit_length'] ** 3
    # data.info["unit_mass"] = data.info["unit_mass"]*fac_len**3
    data.info["unit_time"] = data.info["unit_time"]
    data.info["unit_velocity"] = data.info["unit_velocity"] * fac_len
    data.info["unit_pressure"] = data.info["unit_pressure"] * fac_len ** 2
    data.info["unit_temperature"] = data.info["unit_temperature"] * \
                                    fac_len ** 2 * mu
    units['unit_M_g'] = data.info['unit_mass'].express(C.g)
    units['unit_M_kg'] = data.info['unit_mass'].express(C.kg)
    units['unit_M_Msun'] = data.info['unit_mass'].express(C.Msun)
    units['unit_L_cm'] = data.info['unit_length'].express(C.cm)
    units['unit_L_m'] = data.info['unit_length'].express(C.m)
    units['unit_L_km'] = data.info['unit_length'].express(C.km)
    units['unit_L_au'] = data.info['unit_length'].express(C.au)
    units['unit_L_pc'] = data.info['unit_length'].express(C.pc)
    units['unit_D_cgs'] = data.info['unit_density'].express(C.g / C.cm ** 3)
    units['unit_D_si'] = data.info['unit_density'].express(C.kg / C.m ** 3)
    units['unit_D_Hcc'] = data.info['unit_density'].express(
        C.H_cc) / 0.76 / 2.31
    # ts['unit_D_Hcc'] = data.info['unit_density'].express(C.H_cc)
    units['unit_T_s'] = data.info['unit_time'].express(C.s)
    units['unit_T_yr'] = data.info['unit_time'].express(C.year)
    units['unit_T_Kyr'] = data.info['unit_time'].express(C.year) * (10 ** -3)
    units['unit_T_Myr'] = data.info['unit_time'].express(C.year) * (10 ** -6)
    units['unit_V_cgs'] = data.info['unit_velocity'].express(C.cm / C.s)
    units['unit_V_si'] = data.info['unit_velocity'].express(C.m / C.s)
    units['unit_Pressure'] = data.info['unit_pressure'].express(C.N / C.m ** 2)
    units['unit_Temperature'] = data.info['unit_temperature'].express(C.K)
    print("====== Unit after normalize ======")
    for item in units.items():
        print(item)
    print
    ""
    return units
