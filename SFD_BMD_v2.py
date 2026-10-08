# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 18:04:51 2026

@author: hknar
"""
import numpy as np
import matplotlib.pyplot as plt


def sfd_bmd(weight=3580650, span=61.2, wing_area=437.3, taper=0.22, n=3.5, safety_factor=1.05, num_points=100):

    # ============================================================
    # Aircraft / Wing Data
    # ============================================================
    semi_span = span / 2   # m

    aspect_ratio = span**2 / wing_area

    c_root = 2.0 * wing_area / (span * (1.0 + taper))
    c_tip = taper * c_root

    # ============================================================
    # Design Lift
    # ============================================================
    design_lift = weight * n * safety_factor

    # L/S = q * CL
    q_Cl = design_lift / wing_area

    # ============================================================
    # Spanwise Stations
    # Tip -> Root
    # ============================================================
    y = np.linspace(
        semi_span,
        0.0,
        num_points
    )

    n_stations = len(y)

    # Normalized spanwise coordinate
    eta = 2.0 * y / span

    # ============================================================
    # Trapezoidal Chord Distribution
    # ============================================================
    c_trap = c_root * (
        1.0 - (1.0 - taper) * eta
    )

    # ============================================================
    # Equivalent Elliptical Chord
    # ============================================================
    c_root_ellipse = (
        4.0 * wing_area
        / (np.pi * span)
    )

    c_ellipse = (
        c_root_ellipse
        * np.sqrt(np.maximum(1.0 - eta**2, 0.0))
    )

    # ============================================================
    # Schrenk Average Chord
    # ============================================================
    c_avg = (
        c_trap + c_ellipse
    ) / 2.0

    # ============================================================
    # Local Lift Coefficient Ratio
    # ============================================================
    cl_ratio = c_avg / c_trap

    # ============================================================
    # Initialize Arrays
    # ============================================================
    delta_y = np.zeros(n_stations)
    c_trap_avg = np.zeros(n_stations)
    elem_area = np.zeros(n_stations)
    dy_c_sq = np.zeros(n_stations)
    cl_ratio_avg = np.zeros(n_stations)

    elem_lift = np.zeros(n_stations)
    shear_force = np.zeros(n_stations)
    bending_moment = np.zeros(n_stations)

    distributed_lift = q_Cl * c_avg
    trap = getattr(np, 'trapezoid', np.trapz if hasattr(np, 'trapz') else None)
    distributed_lift *= (design_lift/2) / (-trap(distributed_lift, y))
    for i in range(1, n_stations):
        delta_y[i] = y[i-1] - y[i]
        c_trap_avg[i] = .5*(c_trap[i-1]+c_trap[i])
        elem_area[i] = c_trap_avg[i]*delta_y[i]
        dy_c_sq[i] = delta_y[i]*c_trap_avg[i]**2
        cl_ratio_avg[i] = .5*(cl_ratio[i-1]+cl_ratio[i])
        elem_lift[i] = .5*(distributed_lift[i-1]+distributed_lift[i])*delta_y[i]
        shear_force[i] = shear_force[i-1]+elem_lift[i]
        bending_moment[i] = bending_moment[i-1]+.5*(shear_force[i-1]+shear_force[i])*delta_y[i]

    return {
        "y": y,
        "eta": eta,
        "chord_trapezoidal": c_trap,
        "chord_elliptical": c_ellipse,
        "chord_schrenk": c_avg,
        "cl_ratio": cl_ratio,
        "delta_y": delta_y,
        "element_chord_avg": c_trap_avg,
        "element_area": elem_area,
        "dy_c_sq": dy_c_sq,
        "element_cl_ratio_avg": cl_ratio_avg,
        "element_lift": elem_lift,
        "shear_force": shear_force,
        "bending_moment": bending_moment,
        "design_lift": design_lift,
        "q_Cl": q_Cl,
    }


def plot_distributions(results, save_path=None):

    y = results["y"]

    design_lift = results["design_lift"]
    q_Cl = results["q_Cl"]

    # ============================================================
    # Plot Style
    # ============================================================
    plt.style.use(
        "seaborn-v0_8-whitegrid"
        if "seaborn-v0_8-whitegrid" in plt.style.available
        else "default"
    )

    fig, axs = plt.subplots(
        2,
        2,
        figsize=(14, 10)
    )

    # ============================================================
    # 1. Chord Distributions
    # ============================================================
    axs[0, 0].plot(
        y,
        results["chord_trapezoidal"],
        label="Trapezoidal",
        color="#3498db",
        linewidth=2
    )

    axs[0, 0].plot(
        y,
        results["chord_elliptical"],
        label="Elliptical",
        color="#2ecc71",
        linewidth=2,
        linestyle="--"
    )

    axs[0, 0].plot(
        y,
        results["chord_schrenk"],
        label="Schrenk Average",
        color="#e74c3c",
        linewidth=2.5
    )

    axs[0, 0].set_title(
        "Chord Distribution along Semi-Span",
        fontsize=12,
        fontweight="bold"
    )

    axs[0, 0].set_xlabel(
        "Spanwise position y (m)"
    )

    axs[0, 0].set_ylabel(
        "Chord c (m)"
    )

    axs[0, 0].legend()
    axs[0, 0].grid(
        True,
        linestyle=":",
        alpha=0.6
    )

    # ============================================================
    # 2. Lift Distribution
    # ============================================================
    lift_dist = (
        results["chord_schrenk"]
        * q_Cl
    )

    axs[0, 1].plot(
        y,
        lift_dist,
        color="#9b59b6",
        linewidth=2.5
    )

    axs[0, 1].set_title(
        "Lift Load Distribution L'(y)",
        fontsize=12,
        fontweight="bold"
    )

    axs[0, 1].set_xlabel(
        "Spanwise position y (m)"
    )

    axs[0, 1].set_ylabel(
        "Lift per unit span L' (N/m)"
    )

    axs[0, 1].grid(
        True,
        linestyle=":",
        alpha=0.6
    )

    # ============================================================
    # 3. Shear Force
    # ============================================================
    axs[1, 0].plot(
        y,
        results["shear_force"] / 1e3,
        color="#e67e22",
        linewidth=2.5
    )

    axs[1, 0].set_title(
        "Shear Force Distribution V(y)",
        fontsize=12,
        fontweight="bold"
    )

    axs[1, 0].set_xlabel(
        "Spanwise position y (m)"
    )

    axs[1, 0].set_ylabel(
        "Shear Force (kN)"
    )

    axs[1, 0].grid(
        True,
        linestyle=":",
        alpha=0.6
    )

    # ============================================================
    # 4. Bending Moment
    # ============================================================
    axs[1, 1].plot(
        y,
        results["bending_moment"] / 1e6,
        color="#1abc9c",
        linewidth=2.5
    )

    axs[1, 1].set_title(
        "Bending Moment Distribution M(y)",
        fontsize=12,
        fontweight="bold"
    )

    axs[1, 1].set_xlabel(
        "Spanwise position y (m)"
    )

    axs[1, 1].set_ylabel(
        "Bending Moment (MN·m)"
    )

    axs[1, 1].grid(
        True,
        linestyle=":",
        alpha=0.6
    )

    plt.tight_layout()

    if save_path:
        plt.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight"
        )
        print(f"Plots saved to {save_path}")
    else:
        plt.show()

   # plt.close()


# ================================================================
# Main Program
# ================================================================
if __name__ == "__main__":
    results = sfd_bmd()
    plot_distributions(results)
    print(
        f"Design Lift = "
        f"{results['design_lift'] / 1e6:.3f} MN"
    )

    print(
        f"Root Shear Force = "
        f"{results['shear_force'][-1] / 1e3:.3f} kN"
    )

    print(
        f"Root Bending Moment = "
        f"{results['bending_moment'][-1] / 1e6:.3f} MN.m"
    )

    plot_distributions(results)