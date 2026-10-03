"""
PRDATEN English Summary Report Generator
Replicates Fortran PDRUCK format (PRDATEN.DAT) with standardized English terminology.
"""
import math


def generate_prdaten_report(result: dict, values: dict) -> str:
    tlast = result.get("tlast_result", {})
    pdruck = result.get("pdruck_result", {})
    eexe = pdruck.get("eexe", {})
    nl = pdruck.get("no_load", {})

    rated = tlast["points"][0] if tlast.get("points") else {}

    p_poles = int(float(values.get("p", 2.0)))
    fn = float(values.get("fn", 50.0))
    un = float(values.get("un", 400.0))
    pn = float(values.get("pn", 90.0))
    n_sync = (fn / p_poles) * 60.0

    sf = tlast.get("sfeld", [])
    while len(sf) < 5:
        sf.append(0.0)

    n_rated = (1.0 - sf[3]) * n_sync

    ishalt = int(values.get("ishalt", 1))
    i_scale = 1.0 if ishalt in (0, 2) else math.sqrt(3.0)
    i_rated = rated.get("i1", 0.0) * i_scale

    conn = "DELTA (D)" if ishalt in (1, 3) else "STAR (Y)"

    lines = [
        " " * 7 + "T E S T   B A Y   I N F O R M A T I O N   >      0<",
        "",
        f"       TYPE  : {values.get('reke', ''):<25}  NAME  : {values.get('name', ''):<20}",
        f"       ORDER : {values.get('fbnr', ''):<25}  M-NO. : {values.get('monr', ''):<20}",
        f"       CODE  : {values.get('kennwort', ''):<25}  DATE  : {values.get('datum', ''):<8}",
        "",
        f"       3~ INDUCTION MOTOR  [{conn}]",
        f"         {un:5.0f}. V   {i_rated:7.0f}. A   {pn:7.0f}. kW   {n_rated:7.1f} RPM   {fn:5.1f} HZ",
        "",
        "       P/PN         1/4       2/4       3/4       4/4       5/4",
        ""
    ]

    cf = tlast.get("cfeld", [])
    while len(cf) < 5: cf.append(0.0)
    ef = tlast.get("efeld", [])
    while len(ef) < 5: ef.append(0.0)

    lines.append(f"       COS PHI   {cf[0]:7.3f}   {cf[1]:7.3f}   {cf[2]:7.3f}   {cf[3]:7.3f}   {cf[4]:7.3f}")
    lines.append("")
    lines.append(
        f"       ETA (%)   {ef[0] * 100.:7.2f}   {ef[1] * 100.:7.2f}   {ef[2] * 100.:7.2f}   {ef[3] * 100.:7.2f}   {ef[4] * 100.:7.2f}")
    lines.append("")
    lines.append(f"       SLIP      {sf[0]:7.5f}   {sf[1]:7.5f}   {sf[2]:7.5f}   {sf[3]:7.5f}   {sf[4]:7.5f}")
    lines.append("")
    lines.append("")
    lines.append("       U/UN      IA/IN     MA/MN     MK/MN     STATOR (K/S)  ROTOR (K/S)")
    lines.append("")

    corr = tlast.get("corrected", {})
    dh = corr.get("densities_heating", {})

    ia_in = (corr.get("i1", 0.0) * i_scale) / i_rated if i_rated else 0.0
    ma_mn = corr.get("m", 0.0) / rated.get("m", 1.0) if rated.get("m") else 0.0
    mk_mn = tlast.get("mkipp", 0.0) / rated.get("m", 1.0) if rated.get("m") else 0.0
    lines.append(
        f"       1.00   {ia_in:8.2f}  {ma_mn:8.2f}  {mk_mn:8.2f}    {dh.get('ta1', 0.0):10.2f}    {dh.get('ta2', 0.0):10.2f}")
    lines.append("")

    red = pdruck.get("reduced_voltage")
    if red:
        u_ratio = red.get("u1str_over_u1", 0.0)
        ia_in_r = red.get("i1", 0.0) / i_rated if i_rated else 0.0
        ma_mn_r = red.get("m", 0.0) / rated.get("m", 1.0) if rated.get("m") else 0.0
        mk_mn_r = red.get("kipp_m", 0.0) / rated.get("m", 1.0) if rated.get("m") else 0.0
        lines.append(f"      {u_ratio:5.2f}   {ia_in_r:8.2f}  {ma_mn_r:8.2f}  {mk_mn_r:8.2f}")
        lines.append("")
        lines.append("       (LINE 2 INCLUDES VOLTAGE DROP ACROSS SERIES REACTANCE IF CONFIGURED)")
    lines.append("")
    lines.append("")

    lines.append("       S H O R T - C I R C U I T   C H A R A C T E R I S T I C")
    lines.append("")
    lines.append("            U1(V)     IA(A)    MA(NM)      COSA    PA(KW)   SA(KVA)    KKA")
    for pt in pdruck.get("short_circuit_curve", []):
        lines.append(
            f"          {pt.get('u1', 0):7.1f}   {pt.get('i1', 0):7.1f}   {pt.get('m', 0):7.0f}     {pt.get('cosph', 0):6.3f}   {pt.get('p1_kw', 0):7.1f}  {pt.get('s_kva', 0):8.1f}   {pt.get('kka', 0.84):5.2f}")
    lines.append("")
    lines.append("")

    lines.append("       N O - L O A D   (WARM MACHINE)")
    lines.append("")
    lines.append("            I0(A)      COS0    P0(KW)   PFE(KW)  P_MECH(KW)")
    lines.append(
        f"          {nl.get('i1', 0.):7.2f}     {nl.get('cosph', 0.):6.3f}   {nl.get('p1_kw', 0.):7.3f}   {nl.get('pfe1_kw', 0.):7.3f}   {nl.get('prbg0_kw', 0.):7.3f}")
    lines.append("")
    lines.append("")

    dt1 = tlast.get("dt1", 0.0)
    dt2 = tlast.get("dt2", 0.0)
    thaupt = tlast.get("thaupt", 0.0)
    lines.append(f"       DT1 = {dt1:5.1f} K    DT2 = {dt2:5.1f} K    Residual Field Time Constant: {thaupt:5.2f} S")
    lines.append("")
    r1k = pdruck.get("r1k", 0.0)
    r2k = pdruck.get("r2k", 0.0)
    lines.append(f"       Phase Resistances at 20°C in Ohms: {r1k:7.4f} / {r2k:7.4f}")
    lines.append("")
    holauf = result.get("holauf_result") or {}
    th0 = holauf.get("th", 0.0)
    jmotor = 1.2 * result.get("jaktiv", 0.0)
    lines.append(f"       No-Load Acceleration Time at Rated Voltage: {th0:5.1f} S  at J_motor = {jmotor:8.2f} KG*M2")
    lines.append("")

    if eexe:
        lines.append("       Parameters for Explosion-Proof (EEx-e) Certification:")
        lines.append("             Z1     FW1      M1      N2   Q_BAR(mm2)   KR        a   1+TAU/2")
        lines.append(
            f"          {eexe.get('z1', 0.):5.0f}.  {eexe.get('xsip', 0.):6.3f}   {eexe.get('m1', 3.):5.0f}.  {eexe.get('n2', 0.):5.0f}.   {eexe.get('qstab_mm2', 0.):6.1f}   {eexe.get('k1', 1.):5.3f}  {eexe.get('inv_k42', 0.):6.4f}     {eexe.get('eptauh', 1.):5.3f}")
        lines.append("")

    return "\n".join(lines)
