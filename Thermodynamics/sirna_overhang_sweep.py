#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sirna_overhang_sweep.py

Автоматический прогон ОДНОГО и того же пула кандидатов siRNA по нескольким режимам оверхэнгов:
  - dTdT (энергия считается как UU, т.к. ViennaRNA работает с RNA-алфавитом)
  - UU
  - OMe_UU ("умный" оверхэнг: UU + небольшая ΔΔG/ΔTm поправка)
  - F_UU   ("умный" оверхэнг: UU + небольшая ΔΔG/ΔTm поправка)

Скрипт:
1) читает входной CSV из ваших анализаторов (или любой CSV с колонками sense_sequence/antisense_sequence);
2) пересчитывает ΔG дуплекса для каждого режима оверхэнга (ViennaRNA RNA.duplexfold);
3) (опционально) пересчитывает интегральные scores (guide/bias/passenger/total) по той же логике,
   что и в sirna_mod_scoring_addon.py (без модификаций в теле, только влияние через ΔG_overhang);
4) сохраняет:
   - long-format CSV (по строке на (candidate, overhang_kind))
   - summary CSV (корреляции рангов + топ-сдвиги).

Пример:
  python3 sirna_overhang_sweep.py --in sirna_optimized_results_20260125_120000.csv --outdir ./overhang_sweep

Зависимости:
  - pandas, numpy
  - ViennaRNA python bindings (import RNA)
  - файл sirna_overhangs_ext.py (рядом или в PYTHONPATH)

Запуск:
python3 sirna_overhang_sweep.py \
  --in sirna_optimized_results_YYYYMMDD_HHMMSS.csv \
  --outdir ./overhang_sweep \
  --legacy-orientation
"""

from __future__ import annotations
import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any

import numpy as np
import pandas as pd
import os

from sirna_overhangs_ext import OverhangConfig, energy_duplex_with_overhangs, apply_overhang_corrections


# ----------------------------
# Minimal scoring (from addon)
# ----------------------------

@dataclass
class SweepScoreParams:
    # bands
    tm_good: Tuple[float, float] = (48.0, 60.0)
    tm_ok: Tuple[float, float] = (45.0, 65.0)
    dg_good: Tuple[float, float] = (-35.0, -20.0)
    dg_ok: Tuple[float, float] = (-40.0, -15.0)

    # seed scoring
    seed_tm_opt: float = 21.5
    seed_tm_ok: float = 24.0

    # bias scoring (signed asymmetry)
    bias_opt_range: Tuple[float, float] = (0.5, 4.0)
    bias_ok_range: Tuple[float, float] = (-1.0, 6.0)

    # weights
    w_tm: float = 0.25
    w_dg: float = 0.25
    w_seed: float = 0.35
    w_end: float = 0.15

    w_p_end: float = 0.45
    w_p_mod: float = 0.55  # here we approximate passenger suppression as function of GC end only; keep weight for structure

    # total lambdas
    lambda_guide: float = 0.55
    lambda_bias: float = 0.30
    lambda_pass: float = 0.15


def _band_score(x: float, good: Tuple[float, float], ok: Tuple[float, float]) -> float:
    lo_g, hi_g = good
    lo_o, hi_o = ok
    if lo_g <= x <= hi_g:
        return 100.0
    if lo_o <= x < lo_g:
        return 100.0 * (x - lo_o) / max(1e-9, (lo_g - lo_o))
    if hi_g < x <= hi_o:
        return 100.0 * (hi_o - x) / max(1e-9, (hi_o - hi_g))
    dist = min(abs(x - lo_o), abs(x - hi_o))
    return max(0.0, 40.0 - 5.0 * dist)


def _seed_score(seed_tm: Optional[float], opt: float, ok: float) -> float:
    if seed_tm is None or (isinstance(seed_tm, float) and math.isnan(seed_tm)):
        return 50.0
    if seed_tm <= opt:
        return 100.0
    if seed_tm <= ok:
        return 80.0
    if seed_tm <= ok + 2.0:
        return 60.0
    return 30.0


def _end_pref_score(guide5: str, pass5: str) -> float:
    s = 0.0
    if guide5 in {"A", "U"}:
        s += 50.0
    if pass5 in {"G", "C"}:
        s += 50.0
    return s


def score_bias(asym_signed: Optional[float], p: SweepScoreParams) -> float:
    if asym_signed is None or (isinstance(asym_signed, float) and math.isnan(asym_signed)):
        return 50.0
    a = float(asym_signed)
    opt_lo, opt_hi = p.bias_opt_range
    ok_lo, ok_hi = p.bias_ok_range
    if opt_lo <= a <= opt_hi:
        return 100.0
    if ok_lo <= a < opt_lo:
        return 60.0 + 40.0 * (a - ok_lo) / max(1e-9, (opt_lo - ok_lo))
    if opt_hi < a <= ok_hi:
        return 60.0 + 40.0 * (ok_hi - a) / max(1e-9, (ok_hi - opt_hi))
    dist = min(abs(a - ok_lo), abs(a - ok_hi))
    return max(0.0, 40.0 - 10.0 * dist)


def score_guide(tm: float, dg: float, seed_tm: Optional[float], guide5: str, pass5: str, p: SweepScoreParams) -> float:
    tm_s = _band_score(tm, p.tm_good, p.tm_ok)
    dg_s = _band_score(dg, p.dg_good, p.dg_ok)
    seed_s = _seed_score(seed_tm, p.seed_tm_opt, p.seed_tm_ok)
    end_s = _end_pref_score(guide5, pass5)
    raw = p.w_tm * tm_s + p.w_dg * dg_s + p.w_seed * seed_s + p.w_end * end_s
    return round(min(100.0, max(0.0, raw)), 1)


def score_passenger(pass5: str, p: SweepScoreParams) -> float:
    # lightweight proxy: 5' passenger G/C is good (suppresses passenger loading)
    end = 100.0 if pass5 in {"G", "C"} else 60.0
    return round(end, 1)


def score_total(sg: float, sb: float, sp: float, p: SweepScoreParams) -> float:
    return round(min(100.0, p.lambda_guide * sg + p.lambda_bias * sb + p.lambda_pass * sp), 1)


# ----------------------------
# Sweep logic
# ----------------------------

def _spearman_like_rankcorr(a: pd.Series, b: pd.Series) -> float:
    # Spearman = Pearson on ranks
    ra = a.rank(method="average")
    rb = b.rank(method="average")
    return float(ra.corr(rb))


def main():
    ap = argparse.ArgumentParser(description="Sweep overhang modes and compare ranks")
    ap.add_argument("--in", dest="inp", required=True, help="Input CSV from analyzer (must have sense_sequence/antisense_sequence)")
    ap.add_argument("--outdir", default="./overhang_sweep", help="Output directory")
    ap.add_argument("--kinds", default="dTdT,UU,OMe_UU,F_UU", help="Comma-separated overhang kinds")
    ap.add_argument("--legacy-orientation", action="store_true", help="Use legacy overhang placement (sense: prefix; antisense: suffix). Recommended for comparability.")
    ap.add_argument("--flip-asymmetry", action="store_true", help="Flip sign of thermodynamic_asymmetry before bias scoring")
    ap.add_argument("--no-rescore", action="store_true", help="Only recompute energies; do not compute scores")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.inp)

    for col in ["sense_sequence", "antisense_sequence"]:
        if col not in df.columns:
            raise SystemExit(f"Missing required column: {col}")

    if "tm_celsius" not in df.columns:
        # allow running from a plain sequence list: set tm to NaN
        df["tm_celsius"] = np.nan
    if "duplex_energy_kcal_mol" not in df.columns:
        df["duplex_energy_kcal_mol"] = np.nan
    if "thermodynamic_asymmetry" not in df.columns:
        df["thermodynamic_asymmetry"] = np.nan
    if "seed_tm" not in df.columns:
        df["seed_tm"] = np.nan
    if "pair_id" not in df.columns:
        df["pair_id"] = [f"pair_{i:06d}" for i in range(len(df))]

    kinds = [k.strip() for k in args.kinds.split(",") if k.strip()]
    p = SweepScoreParams()

    long_rows = []
    for _, r in df.iterrows():
        sense = str(r["sense_sequence"]).upper()
        anti = str(r["antisense_sequence"]).upper()
        if not sense or not anti:
            continue

        tm0 = r.get("tm_celsius", np.nan)
        try:
            tm0 = float(tm0)
        except Exception:
            tm0 = np.nan

        asym = r.get("thermodynamic_asymmetry", np.nan)
        if args.flip_asymmetry and not (isinstance(asym, float) and math.isnan(asym)):
            asym = -float(asym)

        seed_tm = r.get("seed_tm", np.nan)
        try:
            seed_tm_val = float(seed_tm)
            if math.isnan(seed_tm_val):
                seed_tm_val = None
        except Exception:
            seed_tm_val = None

        for kind in kinds:
            cfg = OverhangConfig(kind=kind)
            dg = energy_duplex_with_overhangs(sense, anti, cfg, legacy_orientation=args.legacy_orientation)
            meta = apply_overhang_corrections(tm0 if not math.isnan(tm0) else 0.0, dg, cfg)

            out = {
                "pair_id": r["pair_id"],
                "start_pos": r.get("start_pos", np.nan),
                "sense_sequence": sense,
                "antisense_sequence": anti,
                "tm_celsius": tm0,
                "thermodynamic_asymmetry": asym,
                "seed_tm": seed_tm_val,
                "overhang_kind": meta["overhang_kind"],
                "dg_overhang": dg,
            }

            if not args.no_rescore:
                if math.isnan(tm0):
                    # If tm absent, guide score will be computed on a default mid value to keep ranking stable,
                    # but user should prefer supplying tm.
                    tm_for = 55.0
                else:
                    tm_for = tm0 if meta["overhang_kind"] in ("dTdT", "UU") else meta["tm_with_overhang_tm_corr"]

                sg = score_guide(tm_for, dg, seed_tm_val, anti[0], sense[0], p)
                sb = score_bias(asym, p)
                sp = score_passenger(sense[0], p)
                st = score_total(sg, sb, sp, p)

                out.update({
                    "score_guide": sg,
                    "score_bias": sb,
                    "score_passenger": sp,
                    "score_total": st,
                })

            long_rows.append(out)

    long_df = pd.DataFrame(long_rows)
    long_path = os.path.join(args.outdir, "overhang_sweep_long.csv")
    long_df.to_csv(long_path, index=False, encoding="utf-8")

    # Summary: rank correlations vs baseline kind[0]
    summary = []
    baseline_kind = kinds[0]
    if not args.no_rescore and "score_total" in long_df.columns:
        base = long_df[long_df["overhang_kind"] == baseline_kind].copy()
        for kind in kinds:
            cur = long_df[long_df["overhang_kind"] == kind].copy()
            merged = base[["pair_id", "score_total"]].merge(cur[["pair_id", "score_total"]], on="pair_id", suffixes=("_base", "_cur"))
            rho = _spearman_like_rankcorr(merged["score_total_base"], merged["score_total_cur"])
            summary.append({"baseline": baseline_kind, "kind": kind, "spearman_rank_corr": round(rho, 4), "n": int(len(merged))})

        # Also compute top movers between baseline and each kind
        movers_rows = []
        base_rank = base[["pair_id", "score_total"]].copy()
        base_rank["rank_base"] = base_rank["score_total"].rank(ascending=False, method="min")
        for kind in kinds:
            cur = long_df[long_df["overhang_kind"] == kind].copy()
            cur_rank = cur[["pair_id", "score_total"]].copy()
            cur_rank["rank_cur"] = cur_rank["score_total"].rank(ascending=False, method="min")
            m = base_rank.merge(cur_rank, on="pair_id", suffixes=("_base", "_cur"))
            m["rank_delta"] = m["rank_base"] - m["rank_cur"]  # positive = improved rank in cur
            # take top 50 movers
            topm = m.sort_values("rank_delta", ascending=False).head(50)
            for _, rr in topm.iterrows():
                movers_rows.append({
                    "kind": kind,
                    "pair_id": rr["pair_id"],
                    "score_base": rr["score_total_base"],
                    "score_cur": rr["score_total_cur"],
                    "rank_base": int(rr["rank_base"]),
                    "rank_cur": int(rr["rank_cur"]),
                    "rank_delta": int(rr["rank_delta"]),
                })
        movers_df = pd.DataFrame(movers_rows)
        movers_path = os.path.join(args.outdir, "overhang_top_movers.csv")
        movers_df.to_csv(movers_path, index=False, encoding="utf-8")
    else:
        summary.append({"baseline": baseline_kind, "kind": "(scores disabled)", "spearman_rank_corr": np.nan, "n": int(len(long_df))})

    summary_df = pd.DataFrame(summary)
    summary_path = os.path.join(args.outdir, "overhang_sweep_summary.csv")
    summary_df.to_csv(summary_path, index=False, encoding="utf-8")

    print(f"✓ long:    {long_path} ({len(long_df)} rows)")
    print(f"✓ summary: {summary_path} ({len(summary_df)} rows)")
    if not args.no_rescore and "score_total" in long_df.columns:
        print(f"✓ movers:  {os.path.join(args.outdir, 'overhang_top_movers.csv')}")
    print("Done.")


if __name__ == "__main__":
    main()
