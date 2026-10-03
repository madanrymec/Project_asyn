"""
Plotting driver (source: PLOTT 7737-7788, FRAME 5233-5335, HEAD 5379-5420,
GITTER 5356-5378, plus the low-level Calcomp/HP-GL pen-plotter primitives
AXIS/SYMBOL/NUMBER/SCALE/DELTI/LINE/DASHL/DASHP/CIRCL/PLOT/PLOTS/NEWPEN/
NEWPLOT/BMASKE/WHERE/FACTOR, source ~9961-12836 tail) - replaced here with
a single matplotlib chart, NOT a line-by-line port.

*** WHY THIS IS A REWRITE, NOT A PORT ***
Tracing every caller of the plotting system (grep for "CALL PLOTT" and
"CALL FRAME" across the whole file) shows PLOTT is called exactly ONCE
in the entire program, right after SLAST in PROGRAM ASYN4 (source line
779). Despite ~1000 lines of Calcomp pen-plotter machinery (coordinate
transforms, HP-GL pen commands, dashed-line drawing, circle drawing,
axis tick-mark placement via SCALE/DELTI), the whole subsystem exists to
draw exactly ONE chart type: the combined torque-speed and current-speed
characteristic (IACHS=3 in FRAME/HEAD's terms - "Strom-Drehmoment /
Drehzahl-Kennlinie"). Given that, faithfully reimplementing AXIS/SYMBOL/
SCALE/CIRCL/etc as generic low-level primitives (which would only ever
be exercised by this one caller anyway) would be pure busywork with no
practical benefit - matplotlib already does axis scaling, tick
placement, and text rendering better than the 1990s HP-GL emulation did.
This module reproduces what the chart actually SHOWS - reading directly
off FRAME's own axis/curve setup (source 5306-5334) - as a normal,
modern chart instead.

What FRAME actually plots (this is the complete specification, read
directly from the source): X-axis = N/N1 (per-unit speed, 0 to ~1.05).
Left Y-axis = M/MN (per-unit torque): the motor's own torque curve, plus
(if a counter-torque/load curve is configured, MG(0)!=0 - see
winding_functions.mgegen) the counter-torque curve overlaid on the same
axes. Right Y-axis = I1/I1N (per-unit stator current). If a reduced-
voltage/starting-reactor run was also computed (st.mnlred, populated
when KU1 or XNETZ are configured), that curve is overlaid too (source
7780-7787, the CALL LINE calls after CALL FRAME).
"""
from __future__ import annotations

from .state import MachineState
from .winding_functions import mgegen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_torque_speed_curve(reke: str, mn: float, ku1: float, xnetz: float,
                             st: MachineState, output_path: str,
                             title: str | None = None):
    """
    Draws the combined torque-speed / current-speed characteristic
    (the FORTRAN's PLOTT+FRAME with IACHS=3), reading directly from
    st.mnline (populated by load_curves.slast()) and, if a reduced-
    voltage/starting-reactor sweep was run, st.mnlred.

    Returns output_path. Raises ValueError if st.mnline hasn't been
    populated yet (mirrors PLOTT's implicit dependency on SLAST having
    already run - there's no meaningful chart to draw otherwise).
    """
    mnline = st.mnline
    if mnline.nanz == 0:
        raise ValueError(
            "plot_torque_speed_curve: st.mnline is empty - call "
            "load_curves.slast() first (this mirrors PLOTT's dependency "
            "on SLAST having already populated the M-N curve).")

    n_pu, m_pu, i1_pu, mg_pu = [], [], [], []
    for j in range(1, mnline.nanz + 1):
        if mnline.npu[j] < 0.0:
            continue
        n_pu.append(mnline.npu[j])
        m_pu.append(mnline.mpu[j])
        i1_pu.append(mnline.i1pu[j])
        mg_pu.append(mgegen(mnline.npu[j], st) / mn if mn else 0.0)

    order = sorted(range(len(n_pu)), key=lambda k: n_pu[k])
    n_pu = [n_pu[k] for k in order]
    m_pu = [m_pu[k] for k in order]
    i1_pu = [i1_pu[k] for k in order]
    mg_pu = [mg_pu[k] for k in order]

    fig, ax_m = plt.subplots(figsize=(9, 6))
    ax_i = ax_m.twinx()

    ax_m.plot(n_pu, m_pu, color="#2b3a55", linewidth=1.8,
              label="Motor torque M/MN")
    if any(v > 0.0 for v in mg_pu):
        ax_m.plot(n_pu, mg_pu, color="#c0392b", linewidth=1.4,
                  linestyle="--", label="Load torque Mg/MN")
    ax_i.plot(n_pu, i1_pu, color="#1f8a70", linewidth=1.4,
              label="Motor current I1/I1N")

    if ku1 >= 0.001 or xnetz >= 0.000001:
        mnlred = st.mnlred
        if mnlred.nanz > 0:
            n_red = [mnlred.npu[j] for j in range(1, mnlred.nanz + 1)]
            m_red = [mnlred.mpu[j] for j in range(1, mnlred.nanz + 1)]
            i1_red = [mnlred.i1pu[j] for j in range(1, mnlred.nanz + 1)]
            ax_m.plot(n_red, m_red, color="#2b3a55", linewidth=1.2,
                      linestyle=":", label="Motor torque (reduced U1)")
            ax_i.plot(n_red, i1_red, color="#1f8a70", linewidth=1.0,
                      linestyle=":", label="Motor current (reduced U1)")

    ax_m.set_xlabel("N / N1  (per-unit speed)")
    ax_m.set_ylabel("M / MN  (per-unit torque)", color="#2b3a55")
    ax_i.set_ylabel("I1 / I1N  (per-unit current)", color="#1f8a70")
    ax_m.tick_params(axis="y", labelcolor="#2b3a55")
    ax_i.tick_params(axis="y", labelcolor="#1f8a70")
    ax_m.set_xlim(0.0, max(n_pu) * 1.02 if n_pu else 1.05)
    ax_m.set_ylim(bottom=0.0)
    ax_i.set_ylim(bottom=0.0)
    ax_m.grid(True, linestyle="-", linewidth=0.4, alpha=0.5)
    ax_m.axhline(0, color="black", linewidth=0.6)

    chart_title = title or (reke.strip() if reke.strip() else "Machine")
    ax_m.set_title(f"{chart_title}\nTorque / current - speed characteristic",
                    fontsize=12, color="#2b3a55")

    lines_m, labels_m = ax_m.get_legend_handles_labels()
    lines_i, labels_i = ax_i.get_legend_handles_labels()
    ax_m.legend(lines_m + lines_i, labels_m + labels_i, loc="upper left",
                fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    return output_path
