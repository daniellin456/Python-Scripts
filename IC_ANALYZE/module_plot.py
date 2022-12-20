import numpy as np
from matplotlib.colors import *
import matplotlib.pyplot as plt


def plot_equilibrium_temperature(fig, ax, nH_array, balance_temperature):
    ax.plot(nH_array, balance_temperature, linewidth=3, linestyle="--", color="black")
    return


def plot_dLambdadRho_by_np_diff(fig, ax, nH_array, T_array, dLambda_dRho_by_np_diff):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array[1:]), np.log10(dLambda_dRho_by_np_diff[i]), label="T =" + str(T_array[i]))

    ax.set_ylim(-3, 3)
    ax.legend()
    return


def plot_dLambdadRho_by_analytical(fig, ax, nH_array, T_array, dLambda_dRho_by_analytical):
    for i in range(0, len(T_array)):
        ax.plot(np.log10(nH_array), np.log10(dLambda_dRho_by_analytical[i]), label="T =" + str(T_array[i]))

    ax.set_ylim(-3, 3)
    ax.legend()
    return


def plot_dLambdadTemp(fig, ax, nH_array, T_array, dLambda_dTemp):
    specific_Rho_Array = np.asarray([1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e6, 1e6])
    for i in range(0, len(specific_Rho_Array)):
        ax.plot(np.log10(T_array[1:]), dLambda_dTemp[i], label=r"$\rho =$" + str(specific_Rho_Array[i]))
    ax.legend()
    return


def plot_m(fig, ax, extent, m_matrix):
    print("Extrema of m: min:%24.14e, max: %24.14e" % (np.min(m_matrix), np.max(m_matrix)))

    im = ax.imshow(m_matrix.T, origin='lower', cmap="rainbow", norm=LogNorm(vmin=1e-1, vmax=1e2), extent=extent)
    ax.set_title("m(nH, T)")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend='both', label='m')
    return


def plot_n(fig, ax, extent, n_matrix):
    print("Extrema of n: min:%24.14e, max: %24.14e" % (np.min(n_matrix), np.max(n_matrix)))
    im = ax.imshow(n_matrix.T, origin='lower', cmap="rainbow", norm=LogNorm(vmin=1e-2, vmax=1e2), extent=extent)
    ax.set_title("n(nH, T)")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend='both', label='n')
    return


def plot_Gamma1(fig, ax, extent, Gamma1):
    print("Extrema of Gamma1: min:%24.14e, max: %24.14e" % (np.min(Gamma1), np.max(Gamma1)))
    im = ax.imshow(Gamma1.T, origin='lower', cmap="rainbow", norm=SymLogNorm(linthresh=1, vmin=-10, vmax=10),
                   extent=extent)
    ax.set_title(r"$\Gamma_1\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_1$")
    return


def plot_Gamma2(fig, ax, extent, Gamma2):
    print("Extrema of Gamma2: min:%24.14e, max: %24.14e" % (np.min(Gamma2), np.max(Gamma2)))
    im = ax.imshow(Gamma2.T, origin="lower", cmap="rainbow", norm=SymLogNorm(linthresh=1, vmin=-10, vmax=10),
                   extent=extent)
    ax.set_title(r"$\Gamma_2\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_2$")
    return


def plot_Gamma3(fig, ax, extent, Gamma3):
    print("Extrema of Gamma_3: min: %24.14e, max: %24.14e" % (np.min(Gamma3), np.max(Gamma3)))
    im = ax.imshow(Gamma3.T, origin="lower", cmap="rainbow", norm=SymLogNorm(linthresh=0.01, vmin=-10, vmax=10),
                   extent=extent)
    ax.set_title(r"$\Gamma_3\rm{(nH,T)}$")
    ax.set_xlabel(r'$\rm{log}\; nH \; (\rm{cm^{-3}})$')
    ax.set_ylabel(r'$\rm{log}\; T \; (\rm{K})$')
    fig.colorbar(im, ax=ax, extend="both", label=r"$\Gamma_3$")
    return


def plot_Gamma1_m1(fig, ax, nH_mesh, T_mesh, Gamma1_m1):
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma1_m1.T, color="white", density=0.5)
    return


def plot_Gamma2_m1(fig, ax, nH_mesh, T_mesh, Gamma2_m1):
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma2_m1.T, color="white", density=0.5)
    return


def plot_Gamma3_m1(fig, ax, nH_mesh, T_mesh, Gamma3_m1):
    ax.streamplot(nH_mesh, T_mesh, np.ones(shape=nH_mesh.shape), Gamma3_m1.T, color="white", density=0.5)
    return
