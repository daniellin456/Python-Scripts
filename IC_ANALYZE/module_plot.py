from matplotlib.colors import *


def plot_equilibrium_temperature(fig, ax, nH_array, balance_temperature):
    ax.plot(nH_array, balance_temperature, linewidth=3, linestyle="--", color="black")
    return


def plot_dTotal_Lambda_dRho_by_np_diff(fig, ax, nH_array, T_array, dTotal_Lambda_dRho_by_np_diff):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array[1:]), np.log10(dTotal_Lambda_dRho_by_np_diff[i]), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Lambda_{tot}}{d\rho}}\;by\;np.diffl$")
    ax.set_xlabel(r"$\rm{log\;nH\;cm^{-3}}$")
    ax.set_ylabel(r"$\rm{\frac{d\Lambda_{tot}}{d\rho}}$")
    ax.set_ylim(-3, 3)
    ax.legend()
    return


def plot_dTotal_Lambda_dRho_by_analytical(fig, ax, nH_array, T_array, dTotal_Lambda_dRho_by_analytical):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array), np.log10(dTotal_Lambda_dRho_by_analytical[i]), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Lambda_{tot}}{d\rho}}\;by\;analytical$")
    ax.set_xlabel(r"$\rm{log\;nH\;cm^{-3}}$")
    ax.set_ylabel(r"$\rm{\frac{d\Lambda_{tot}}{d\rho}}$")
    ax.set_ylim(-3, 3)
    ax.legend()
    return


def plot_m(fig, ax, extent, m_matrix):
    print("Extrema of m: min:%24.14e, max: %24.14e" % (np.min(m_matrix), np.max(m_matrix)))

    im = ax.imshow(m_matrix.T, origin='lower', cmap="rainbow", vmin=0.1, vmax=3, extent=extent)
    ax.set_title("m(nH, T)")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend='both', label='m')
    return


def plot_m_contour(fig, ax, nH_mesh, T_mesh, m_matrix):
    CS = ax.contour(nH_mesh, T_mesh, m_matrix.T, levels=[2])
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    return


def plot_n(fig, ax, extent, n_matrix):
    print("Extrema of n: min:%24.14e, max: %24.14e" % (np.min(n_matrix), np.max(n_matrix)))
    im = ax.imshow(n_matrix.T, origin='lower', cmap="rainbow", vmin=0.1, vmax=2, extent=extent)
    ax.set_title("n(nH, T)")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend='both', label='n')
    return


def plot_n_contour(fig, ax, nH_mesh, T_mesh, n_matrix):
    CS = ax.contour(nH_mesh, T_mesh, n_matrix.T, levels=[0.5, 1])
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    return


def plot_Gamma1(fig, ax, extent, Gamma1):
    print("Extrema of Gamma1: min:%24.14e, max: %24.14e" % (np.min(Gamma1), np.max(Gamma1)))
    im = ax.imshow(Gamma1.T, origin='lower', cmap="rainbow", vmin=-3, vmax=3, extent=extent)
    ax.set_title(r"$\Gamma_1\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_1$", ticks=np.arange(-3, 3.5, 0.5))
    return


def plot_Gamma2(fig, ax, extent, Gamma2):
    print("Extrema of Gamma2: min:%24.14e, max: %24.14e" % (np.min(Gamma2), np.max(Gamma2)))
    im = ax.imshow(Gamma2.T, origin="lower", cmap="rainbow", vmin=-3, vmax=3, extent=extent)
    ax.set_title(r"$\Gamma_2\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_2$", ticks=np.arange(-3, 3.5, 0.5))
    return


def plot_Gamma3(fig, ax, extent, Gamma3):
    print("Extrema of Gamma3: min: %24.14e, max: %24.14e" % (np.min(Gamma3), np.max(Gamma3)))
    im = ax.imshow(Gamma3.T, origin="lower", cmap="rainbow", vmin=-3, vmax=3, extent=extent)
    ax.set_title(r"$\Gamma_3\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_3$", ticks=np.arange(-3, 3.5, 0.5))
    return


def plot_Gamma1_m1(fig, ax, nH_mesh, T_mesh, Gamma1_m1):
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma1_m1.T, color="white", density=1, arrowsize=0)
    return


def plot_Gamma2_m1(fig, ax, nH_mesh, T_mesh, Gamma2_m1):
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma2_m1.T, color="white", density=1, arrowsize=0)
    return


def plot_Gamma3_m1(fig, ax, nH_mesh, T_mesh, Gamma3_m1):
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma3_m1.T, color="white", density=1, arrowsize=0)
    return


def plot_Gamma2_contour(fig, ax, nH_mesh, T_mesh, Gamma2):
    CS = ax.contour(nH_mesh, T_mesh, Gamma2.T, levels=[1, 4 / 3, 7 / 5, 5 / 3, 1.9, 2.0])
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    return


def plot_dLambda_dRho_by_np_diff(fig, ax, nH_array, T_array, dLambda_dRho_by_np_diff):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array[1:]), np.log10(np.abs(dLambda_dRho_by_np_diff[i])), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Lambda}{d\rho}}\;by\;np.diff$")
    ax.set_xlabel(r"$\rm{log\;nH\;cm^{-3}}$")
    ax.set_ylabel(r"$\rm{\frac{d\Lambda}{d\rho}}$")
    ax.legend()
    return


def plot_dLambda_dRho_by_analytical(fig, ax, nH_array, T_array, dLambda_dRho_by_analytical):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array), np.log10(np.abs(dLambda_dRho_by_analytical[i])), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Lambda}{d\rho}}\;by\;analytical$")
    ax.set_xlabel(r"$\rm{log\;nH\;cm^{-3}}$")
    ax.set_ylabel(r"$\rm{\frac{d\Lambda}{d\rho}}$")
    ax.legend()
    return


def plot_dGamma_dRho_by_np_diff(fig, ax, nH_array, T_array, dGamma_dRho_by_np_diff):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array[1:]), np.log10(dGamma_dRho_by_np_diff[i]), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Gamma}{d\rho}}\;by\;np.diff$")
    ax.set_xlabel(r"$\rm{log\;nH\;cm^{-3}}$")
    ax.set_ylabel(r"$\rm{\frac{d\Gamma}{d\rho}}$")
    ax.set_ylim(-36, -24)

    ax.legend()
    return


def plot_dGamma_dRho_by_analytical(fig, ax, nH_array, T_array, dLambda_dRho_by_analytical):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array), np.log10(dLambda_dRho_by_analytical[i]), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Gamma}{d\rho}}\;by\;analytical$")
    ax.set_xlabel(r"$\rm{log\;nH\;cm^{-3}}$")
    ax.set_ylabel(r"$\rm{\frac{d\Gamma}{d\rho}}$")
    ax.set_ylim(-36, -24)

    ax.legend()
    return


def plot_ionization(fig, ax, nH_array, T_array, ionization):
    ax_ionization = ax.twinx()
    for i in range(0, len(T_array)):
        ax_ionization.plot(np.log10(nH_array), np.log10(ionization[i]), label="T =" + str(T_array[i]), linestyle='--')
    ax_ionization.set_ylabel('log10 ionization')
    ax_ionization.legend(loc="center right")
    return


def plot_individual_cooling_heating(ax, nH_array, T, CII_cooling_rate, OI_cooling_rate, H2_cooling_rate,
                                    CII_line_cooling_rate, OI_line_cooling_rate, Recombination_cooling_rate,
                                    UV_heating_rate):
    ax.plot(np.log10(nH_array), np.log10(CII_cooling_rate), label="CII", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(OI_cooling_rate), label="OI", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(H2_cooling_rate), label="H2", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(CII_line_cooling_rate), label="CII line", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(OI_line_cooling_rate), label="OI Line", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(Recombination_cooling_rate), label="Rec", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(UV_heating_rate), label="UV")
    ax.set_title("T = " + str(T) + "K")
    ax.set_xlabel(r"$\rm{log10 \; nH \; (cm^{-3})}$")
    ax.set_ylabel(r"$\rm{log10} \; n\Gamma, \rm{log10} \; n^2 \Lambda$")
    ax.legend(loc='lower right')
    return
