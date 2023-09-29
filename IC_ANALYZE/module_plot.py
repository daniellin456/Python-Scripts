from matplotlib.colors import *
from matplotlib.ticker import MultipleLocator


def plot_equilibrium_temperature(fig, ax, nH_array, balance_temperature):
    ax.plot(nH_array, balance_temperature, linewidth=2, linestyle="--", color="black")
    return


def plot_total_cooling_rate(fig, ax, extent, total_cooling_power):
    print("Extrema of total_cooling_power: min:%24.14e, max: %24.14e" % (
        np.min(total_cooling_power), np.max(total_cooling_power)))
    im = ax.imshow(total_cooling_power.T, origin="lower", cmap="bwr", extent=extent, interpolation="none",
                   norm=SymLogNorm(linthresh=1e-26, vmin=-1e-23, vmax=1e-23))
    ax.set_title("Cooling and Heating")
    ax.set_xlabel(r'$\rm{log_{10}}\; n_H \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log_{10}}\; T \; (\rm{K})$')
    ax.set_xticks([-2, 0, 2, 4, 6, 8, 10])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_yticks([0, 1, 2, 3, 4])
    ticks = [-1e-23, -1e-24, -1e-25, -1e-26, 0, 1e-26, 1e-25, 1e-24, 1e-23]
    cb = fig.colorbar(im, ax=ax, extend='both', label=r"$\rm{\Gamma - n_H\Lambda} \;(erg \;s^{-1})$", location="bottom",
                      pad=0.3, ticks=ticks)
    cb.ax.tick_params(labelsize=10)
    cb.ax.minorticks_on()
    return


def plot_dTotal_Lambda_dRho_by_np_diff(fig, ax, nH_array, T_array, dTotal_Lambda_dRho_by_np_diff):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array[1:]), np.log10(dTotal_Lambda_dRho_by_np_diff[i]), label="T =" + str(T_array[i]))

    ax.set_title(r"$\rm{\frac{d\Lambda_{tot}}{d\rho}}\;by\;np.diff$")
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

    im = ax.imshow(m_matrix.T, origin='lower', cmap="rainbow", vmin=0.1, vmax=3, extent=extent, interpolation="none")
    ax.set_title(r"$\rm{m(n_H, T)}$")
    ax.set_xlabel(r'$\rm{log_{10}}\; n_H \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log_{10}}\; T \; (\rm{K})$')
    ax.set_xticks([-2, 0, 2, 4, 6, 8, 10])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_yticks([0, 1, 2, 3, 4])
    ticks = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    cb = fig.colorbar(im, ax=ax, extend='both', label='m', location="bottom", pad=0.3, ticks=ticks)
    cb.ax.minorticks_on()
    return


def plot_m_contour(fig, ax, nH_mesh, T_mesh, m_matrix):
    CS = ax.contour(nH_mesh, T_mesh, m_matrix.T, levels=[1.5, 2], colors=('grey', 'brown'))
    ax.clabel(CS, CS.levels, inline=True, fontsize=14)
    return


def plot_n(fig, ax, extent, n_matrix):
    print("Extrema of n: min:%24.14e, max: %24.14e" % (np.min(n_matrix), np.max(n_matrix)))
    im = ax.imshow(n_matrix.T, origin='lower', cmap="rainbow", vmin=0.1, vmax=2, extent=extent, interpolation="none")
    ax.set_title(r"$\rm{n(n_H, T)}$")
    ax.set_xlabel(r'$\rm{log_{10}}\; n_H \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log_{10}}\; T \; (\rm{K})$')
    ax.set_xticks([-2, 0, 2, 4, 6, 8, 10])
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.set_yticks([0, 1, 2, 3, 4])
    ticks = [0.1, 0.5, 1.0, 1.50, 2.0]
    cb = fig.colorbar(im, ax=ax, extend='both', label='n', location="bottom", pad=0.3, ticks=ticks)
    cb.ax.minorticks_on()
    return


def plot_n_contour(fig, ax, nH_mesh, T_mesh, n_matrix):
    CS = ax.contour(nH_mesh, T_mesh, n_matrix.T, levels=[0.5, 1], colors=('grey', 'yellow'))
    ax.clabel(CS, CS.levels, inline=True, fontsize=14)
    return


def plot_Gamma1(fig, ax, extent, Gamma1):
    print("Extrema of Gamma1: min:%24.14e, max: %24.14e" % (np.min(Gamma1), np.max(Gamma1)))
    im = ax.imshow(Gamma1.T, origin='lower', cmap="bwr", vmin=-3, vmax=3, extent=extent, interpolation="none")
    ax.set_title(r"$\Gamma_1\rm{(n_H,T)}$")
    ax.set_xlabel(r'$\rm{log}\; n_H \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    cb = fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_1$", ticks=np.arange(-3, 3.5, 0.5), location="bottom",
                      pad=0.3)
    cb.ax.tick_params(labelsize=14)
    return


def plot_Gamma2(fig, ax, extent, Gamma2):
    print("Extrema of Gamma2: min:%24.14e, max: %24.14e" % (np.min(Gamma2), np.max(Gamma2)))
    im = ax.imshow(Gamma2.T, origin="lower", cmap="bwr", extent=extent, interpolation="none",
                   norm=TwoSlopeNorm(vcenter=1.0, vmin=-3, vmax=3))
    # ax.set_title(r"$\Gamma \; \rm{(nH,T)}, \;\; \Gamma\rm{=1+(3-2m)/(2n-2)}$")
    ax.set_title(r"$\gamma\rm{(n_H,T)}$")
    ax.set_xlabel(r'$\rm{log_{10}}\; n_H \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log_{10}}\; T \; (\rm{K})$')
    ax.set_xticks([-2, 0, 2, 4, 6, 8, 10])
    ax.set_yticks([0, 1, 2, 3, 4])
    # fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma$", ticks=np.arange(-3, 3.5, 0.5), location="bottom")
    cb = fig.colorbar(im, ax=ax, extend="both", label=r"$\gamma$", location="bottom", pad=0.3)
    cb.ax.tick_params(labelsize=14)
    cb.ax.minorticks_on()
    return


def plot_Gamma3(fig, ax, extent, Gamma3):
    print("Extrema of Gamma3: min: %24.14e, max: %24.14e" % (np.min(Gamma3), np.max(Gamma3)))
    im = ax.imshow(Gamma3.T, origin="lower", cmap="bwr", vmin=-3, vmax=3, extent=extent, interpolation="none")
    ax.set_title(r"$\Gamma_3\rm{(n_H,T)}$")
    ax.set_xlabel(r'$\rm{log}\; n_H \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    cb = fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_3$", ticks=np.arange(-3, 3.5, 0.5), location="bottom",
                      pad=0.3)
    cb.ax.tick_params(labelsize=8)
    return im


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
    CS = ax.contour(nH_mesh, T_mesh, Gamma2.T, levels=[1.67])
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


def plot_individual_cooling_rate_vs_nH(ax, nH_array, T, CII_cooling_rate, OI_cooling_rate, H2_cooling_rate,
                                       CII_line_cooling_rate, OI_line_cooling_rate,
                                       Recombination_cooling_rate,
                                       UV_heating_rate):
    ax.plot(np.log10(nH_array), np.log10(CII_cooling_rate), label="CII", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(OI_cooling_rate), label="OI", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(H2_cooling_rate), label=r"$\rm{H_2}$", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(CII_line_cooling_rate), label="CII metastable Line", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(OI_line_cooling_rate), label="OI metastable Line", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(Recombination_cooling_rate), label="Rec", linestyle="--")
    ax.plot(np.log10(nH_array), np.log10(UV_heating_rate), label="UV")
    ax.set_ylim(-30, 0)
    ax.set_title("T = " + str(T) + "K", fontsize=16)
    # ax.set_xlabel(r"$\rm{log_{10} \; n_H \; (cm^{-3})}$", fontsize=10)
    # ax.set_ylabel(r"$\rm{log_{10}} \; n\Gamma, \rm{log_{10}} \; n^2 \Lambda \; (\rm{ergs \; cm^{-3}\; s^{-1}})$",
    #              fontsize=10)
    ax.legend(loc="upper left", fontsize=10)
    return


def plot_individual_cooling_rate_vs_T(ax, nH, T_array, CII_cooling_rate, OI_cooling_rate, H2_cooling_rate,
                                      CII_line_cooling_rate, OI_line_cooling_rate, Recombination_cooling_rate,
                                      UV_heating_rate):
    ax.plot(np.log10(T_array), np.log10(CII_cooling_rate), label="CII", linestyle="--")
    ax.plot(np.log10(T_array), np.log10(OI_cooling_rate), label="OI", linestyle="--")
    ax.plot(np.log10(T_array), np.log10(H2_cooling_rate), label=r"$\rm{H_2}$", linestyle="--")
    ax.plot(np.log10(T_array), np.log10(CII_line_cooling_rate), label="CII metastable Line", linestyle="--")
    ax.plot(np.log10(T_array), np.log10(OI_line_cooling_rate), label="OI metastable Line", linestyle="--")
    ax.plot(np.log10(T_array), np.log10(Recombination_cooling_rate), label="Rec", linestyle="--")
    ax.plot(np.log10(T_array), np.log10(UV_heating_rate), label="UV", linestyle="-")
    ax.set_title(r"$\rm{n_H} = $" + "{:.2e}".format(nH) + r" $\rm{cm^{-3}}$", fontsize=16)
    ax.set_ylim(-30, 0)
    # ax.set_xlabel(r"$\rm{log_{10} \; T \; (K)}$", fontsize=10)
    # ax.set_ylabel(r"$\rm{log_{10}} \; n\Gamma, \rm{log_{10}} \; n^2 \Lambda \; (\rm{ergs \; cm^{-3}\; s^{-1}})$",
    #              fontsize=10)
    ax.legend(loc="upper left", fontsize=10)
    return


def plot_slope_line(ax, slope, intercepts, color):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_vals = np.array(xlim)
    for intercept in intercepts:
        y_vals = intercept + slope * x_vals
        ax.plot(x_vals, y_vals, linestyle=':', color=color)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return


def plot_ax3(fig, ax, extent, data):
    print("Extrema of 3-2a: min: %24.14e, max: %24.14e" % (np.min(data), np.max(data)))
    im = ax.imshow(data.T, origin="lower", cmap="bwr", vmin=-3, vmax=3, extent=extent, interpolation="none")
    ax.set_title(r"$\rm{3-2a}\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\rm{3-2a}$", ticks=np.arange(-3, 3.5, 0.5), location="bottom")
    return


def plot_ax3_contour(fig, ax, nH_mesh, T_mesh, data):
    CS = ax.contour(nH_mesh, T_mesh, data.T, levels=[0])
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    return


def plot_ax4(fig, ax, extent, data):
    print("Extrema of c/(\gamma-1)-d: min: %24.14e, max: %24.14e" % (np.min(data), np.max(data)))
    im = ax.imshow(data.T, origin="lower", cmap="bwr", vmin=-3, vmax=3, extent=extent, interpolation="none")
    ax.set_title(r"$\rm{c/(gamma-1)-d}\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$c/(\gamma-1)-d$", ticks=np.arange(-3, 3.5, 0.5), location="bottom")
    return


def plot_ax4_contour(fig, ax, nH_mesh, T_mesh, data):
    CS = ax.contour(nH_mesh, T_mesh, data.T, levels=[0])
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    return


def plot_ax5(fig, ax, extent, data):
    print("Extrema of  c/(gamma-1)+b+2: min: %24.14e, max: %24.14e" % (np.min(data), np.max(data)))
    im = ax.imshow(data.T, origin="lower", cmap="bwr", vmin=-3, vmax=3, extent=extent, interpolation="none")
    ax.set_title(r"$\rm{c/(\gamma-1)+b+2}\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$c/(\gamma-1)+b+2$", ticks=np.arange(-3, 3.5, 0.5),
                 location="bottom")
    return


def plot_ax5_contour(fig, ax, nH_mesh, T_mesh, data):
    CS = ax.contour(nH_mesh, T_mesh, data.T, levels=[0])
    ax.clabel(CS, CS.levels, inline=True, fontsize=10)
    return
