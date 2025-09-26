
import streamlit as st
import pandas as pd
import itertools
import math
from io import StringIO

PROTON_MASS = 1.007276466812

st.set_page_config(page_title="Antibody–Peptide Mass Matcher", page_icon="🧪", layout="wide")

st.title("🧪 Antibody–Peptide Mass Matcher")
st.write("Match observed masses (neutral or m/z) to predicted antibody–peptide conjugation stoichiometries.")

with st.sidebar:
    st.header("Inputs")
    mab_mass = st.number_input("Parental antibody mass (Da)", min_value=0.0, value=149650.0, step=10.0, format="%.6f")
    delta_per_conj = st.number_input("Δ mass per conjugation (Da)", value=0.0, step=0.001, format="%.6f",
                                     help="Use e.g. -18.01056 for water loss per coupling")
    max_total = st.number_input("Max total copies across all peptides", min_value=1, max_value=64, value=8, step=1)
    charges_txt = st.text_input("Charge states for predicted m/z (optional)", value="",
                                help="e.g., 2 3 4. Leave blank to skip m/z predictions")

    tol_mode = st.radio("Tolerance mode", ["ppm", "Da"], horizontal=True, index=0)
    tol_value = st.number_input("Tolerance value", min_value=0.0, value=20.0 if tol_mode=="ppm" else 1.0,
                                step=0.1, format="%.6f")
    st.caption("Tolerance applies to the space of the values you provide in 'Observed' (neutral masses or m/z).")

    st.markdown("---")
    st.subheader("Upload CSVs")
    st.caption("**Peptides CSV** needs columns: `name`, `peptide_mass` (required). Optional: `linker_mass`, `max_copies`.")
    peptides_file = st.file_uploader("Upload peptides.csv", type=["csv"], key="peps")
    st.caption("**Observed CSV** needs column: `observed_mass` (neutral mass or m/z).")
    observed_file = st.file_uploader("Upload observed.csv", type=["csv"], key="obs")
    run_btn = st.button("Run matching")

def read_peptides(df: pd.DataFrame) -> pd.DataFrame:
    lower_cols = {c.lower(): c for c in df.columns}
    required = {"name", "peptide_mass"}
    if not required.issubset(set(map(str.lower, df.columns))):
        st.error("Peptides CSV must include columns: name, peptide_mass (case-insensitive).")
        return None
    out = df.rename(columns={lower_cols.get("name","name"): "name",
                             lower_cols.get("peptide_mass","peptide_mass"): "peptide_mass"})
    if "linker_mass" in lower_cols:
        out = out.rename(columns={lower_cols["linker_mass"]: "linker_mass"})
    else:
        out["linker_mass"] = 0.0
    if "max_copies" in lower_cols:
        out = out.rename(columns={lower_cols["max_copies"]: "max_copies"})
    else:
        out["max_copies"] = None
    # Ensure numeric
    out["peptide_mass"] = pd.to_numeric(out["peptide_mass"], errors="coerce")
    out["linker_mass"] = pd.to_numeric(out["linker_mass"], errors="coerce").fillna(0.0)
    return out[["name","peptide_mass","linker_mass","max_copies"]]

def generate_combinations(peps: pd.DataFrame, max_total: int):
    ranges = []
    for _, row in peps.iterrows():
        cap = int(row["max_copies"]) if pd.notna(row["max_copies"]) else max_total
        ranges.append(range(0, cap+1))
    combos = []
    for combo in itertools.product(*ranges):
        if sum(combo) == 0:
            continue
        if sum(combo) <= max_total:
            combos.append(combo)
    return combos

def combo_mass(peps: pd.DataFrame, combo, delta_per_conj: float) -> float:
    total = 0.0
    for i, n in enumerate(combo):
        if n == 0:
            continue
        row = peps.iloc[i]
        unit = float(row["peptide_mass"]) + float(row["linker_mass"]) + delta_per_conj
        total += n * unit
    return total

def to_mz(neutral_mass: float, z: int) -> float:
    return (neutral_mass + z * PROTON_MASS) / z

def ppm_error(obs: float, calc: float) -> float:
    return 1e6 * (obs - calc) / calc

def build_predicted(peps: pd.DataFrame, mab_mass: float, delta_per_conj: float, max_total: int, charges):
    combos = generate_combinations(peps, max_total)
    predicted = []
    for combo in combos:
        add_mass = combo_mass(peps, combo, delta_per_conj)
        neutral = mab_mass + add_mass
        stoich = {f"{peps.iloc[i]['name']}_copies": int(n) for i, n in enumerate(combo)}
        stoich["total_copies"] = int(sum(combo))
        stoich["mode"] = "neutral"
        predicted.append((neutral, stoich.copy()))
        for z in charges:
            if z > 0:
                mz = to_mz(neutral, z)
                stoich_z = stoich.copy()
                stoich_z["mode"] = f"m/z (z={z})"
                predicted.append((mz, stoich_z))
    return predicted

def match_masses(observed, predicted, tol_value, tol_mode):
    rows = []
    for obs in observed:
        for mass, stoich in predicted:
            if tol_mode == "Da":
                ok = abs(obs - mass) <= tol_value
            else:
                ok = abs(ppm_error(obs, mass)) <= tol_value
            if ok:
                rows.append({
                    "observed": obs,
                    "predicted": mass,
                    "error_Da": obs - mass,
                    "error_ppm": ppm_error(obs, mass),
                    **stoich
                })
    return pd.DataFrame(rows) if rows else pd.DataFrame()

def parse_charges(ch_txt: str):
    vals = []
    for tok in ch_txt.strip().split():
        try:
            vals.append(int(tok))
        except:
            pass
    return vals

col1, col2 = st.columns([1,1])

with col1:
    if peptides_file is None:
        st.info("Upload **peptides.csv** to begin. You can use the example in the right column.")
    else:
        try:
            pep_df_raw = pd.read_csv(peptides_file)
            peps = read_peptides(pep_df_raw)
            if peps is not None:
                st.subheader("Peptides")
                st.dataframe(peps, use_container_width=True)
        except Exception as e:
            st.error(f"Failed to read peptides CSV: {e}")
            peps = None

with col2:
    if observed_file is None:
        st.info("Upload **observed.csv** to begin. You can use the example below.")
    else:
        try:
            obs_df = pd.read_csv(observed_file)
            if "observed_mass" not in map(str.lower, obs_df.columns):
                st.error("Observed CSV must have a column named 'observed_mass'.")
                obs_df = None
            else:
                # Normalize column name
                real_col = [c for c in obs_df.columns if c.lower()=="observed_mass"][0]
                obs_df = obs_df.rename(columns={real_col:"observed_mass"})
                obs_df["observed_mass"] = pd.to_numeric(obs_df["observed_mass"], errors="coerce")
                st.subheader("Observed")
                st.dataframe(obs_df, use_container_width=True)
        except Exception as e:
            st.error(f"Failed to read observed CSV: {e}")
            obs_df = None

st.markdown("---")

if run_btn:
    if peps is None or obs_df is None:
        st.warning("Please upload both peptides and observed CSVs.")
    else:
        charges = parse_charges(charges_txt)
        predicted = build_predicted(peps, mab_mass, delta_per_conj, int(max_total), charges)
        observed_vals = obs_df["observed_mass"].dropna().astype(float).tolist()
        res = match_masses(observed_vals, predicted, tol_value, tol_mode)
        if res.empty:
            st.warning("No matches found within tolerance.")
        else:
            # Order by observed, then mode, then abs error_ppm (if available)
            if "error_ppm" in res.columns:
                res = res.sort_values(by=["observed","mode","error_ppm"], key=lambda s: s.abs() if s.name=="error_ppm" else s)
            st.success(f"Matched {len(res)} rows.")
            st.dataframe(res, use_container_width=True)

            # Offer download
            csv_bytes = res.to_csv(index=False).encode("utf-8")
            st.download_button("Download results CSV", data=csv_bytes, file_name="mass_matches.csv", mime="text/csv")

# Example section at the bottom
st.markdown("---")
st.subheader("Example CSVs")
st.write("Use these headers and sample rows to prepare your files.")

st.code("""name,peptide_mass,linker_mass,max_copies
PepA,1998.987,1200.123,8
PepB,2999.020,0,4
""", language="csv")

st.code("""observed_mass
155800.1
158900.0
""", language="csv")
