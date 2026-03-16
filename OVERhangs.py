import pandas as pd
import itertools
import RNA

df = pd.read_csv("sirna_merged_all.csv")

OVERHANGS = ["none", "UU", "dTdT", "AA", "GC"]
MODIFICATIONS = ["none", "2OMe", "LNA", "F"]

LIT_CORRECTIONS = {
    "none": (0.0, 0.0),
    "UU": (-1.5, +5.0),
    "dTdT": (-1.2, +4.0),
    "AA": (-0.5, +2.0),
    "GC": (-1.0, +3.0),
    "2OMe": (-0.8, +3.0),
    "LNA": (-2.0, +6.0),
    "F": (-1.3, +4.5)
}

rows = []
for idx, row in df.iterrows():
    seq = row["sense_sequence"]
    name = row["fragment_id"]

    for overhang, mod in itertools.product(OVERHANGS, MODIFICATIONS):
        full_seq = seq + (overhang.replace("dT","T") if overhang != "none" else "")
        fc = RNA.fold_compound(full_seq)
        dG_seq = fc.mfe()[1]

        dG_corr, Tm_corr = LIT_CORRECTIONS.get(overhang, (0,0))
        dG_mod_corr, Tm_mod_corr = LIT_CORRECTIONS.get(mod, (0,0))
        total_dG = dG_seq + dG_corr + dG_mod_corr
        total_Tm = 37 + total_dG + Tm_corr + Tm_mod_corr

        rows.append({
            "Name": name,
            "Overhang": overhang,
            "Modification": mod,
            "Sequence_with_overhang": full_seq,
            "dG_total": round(total_dG,2),
            "Tm_total": round(total_Tm,2)
        })

df_results = pd.DataFrame(rows)
df_results.to_csv("OVERhangs.csv", index=False)
print(df_results.head())