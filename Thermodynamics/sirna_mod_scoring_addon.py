
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sirna_mod_scoring_addon.py

Дополнительный скрипт к вашим siRNA_Thermo_Analyzer*.py.
Задачи:
1) Формально корректная (позиционно-взвешенная) модель модификаций 2'-F vs 2'-O-Me.
2) Раздельный scoring для guide (antisense) и passenger (sense) + bias (strand selection).
3) Экспорт расширенного CSV с полями:
   - tm_corr_total / dg_corr_total
   - tm_corr_guide / dg_corr_guide, tm_corr_passenger / dg_corr_passenger (если включено)
   - score_guide, score_passenger, score_bias, score_total
   - mod_pattern_guide, mod_pattern_passenger (JSON-строки)

Ожидаемый входной CSV:
- из ваших анализаторов, содержащий как минимум колонки:
  sense_sequence, antisense_sequence,
  tm_celsius, duplex_energy_kcal_mol,
  thermodynamic_asymmetry (желательно signed, как во 2-й версии),
  seed_tm (желательно),
  start_pos, pair_id (опционально)

Пример запуска:
  python3 sirna_mod_scoring_addon.py --in sirna_optimized_results_YYYYMMDD_HHMMSS.csv --out sirna_scored_addon.csv

Примечание:
- Скрипт НЕ требует ViennaRNA для расчётов (использует уже посчитанные tm/energy из CSV).
- Все параметры и паттерны модификаций настраиваются в JSON (см. --params-json).
"""

from __future__ import annotations
import argparse
import json
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional, Any

import pandas as pd
import numpy as np


# ----------------------------
# Utility: zones and weights
# ----------------------------

def _zones(L: int) -> Dict[str, List[int]]:
    """Возвращает 1-based индексы позиций по функциональным зонам guide/passenger."""
    seed = list(range(2, min(8, L) + 1))            # 2..8
    center = [i for i in range(9, min(12, L) + 1)]  # 9..12 (включая 10-11 область)
    guide_3p = list(range(max(13, 1), L + 1))       # 13..L
    ends_5p = [1]
    ends_3p = [L]
    return {
        "seed": seed,
        "center": center,
        "guide_3p": guide_3p,
        "end5": ends_5p,
        "end3": ends_3p,
    }


def _default_zone_weights() -> Dict[str, Dict[str, float]]:
    """
    Базовые веса зон (можно калибровать).
    strand: 'guide' (antisense) или 'passenger' (sense)
    zone: seed/center/guide_3p/end5/end3/other
    """
    return {
        "guide": {
            "seed": 1.20,
            "center": 0.70,
            "guide_3p": 1.10,
            "end5": 1.00,
            "end3": 1.10,
            "other": 1.00,
        },
        "passenger": {
            "seed": 1.00,
            "center": 1.20,
            "guide_3p": 1.10,  # для passenger это просто 3'-часть, оставляем 1.1
            "end5": 1.20,
            "end3": 1.20,
            "other": 1.10,
        }
    }


def _zone_of_pos(pos1: int, L: int) -> str:
    z = _zones(L)
    for zone, positions in z.items():
        if pos1 in positions:
            return zone
    return "other"


# ----------------------------
# Modification pattern model
# ----------------------------

@dataclass
class ModParams:
    # Per-modification base effects (before zone weighting)
    # Units: alpha -> °C, beta -> kcal/mol (negative stabilizes)
    alpha_F: float = 1.0
    alpha_OMe: float = 0.7
    beta_F: float = -0.5
    beta_OMe: float = -0.3

    # Zone weights
    zone_weights: Dict[str, Dict[str, float]] = None

    # Optional: apply separate corr for guide/passenger or only total
    separate_corr: bool = True

    # Penalty for central over-modification on guide (positions 9-12)
    center_penalty_per_mod: float = 4.0  # points per mod in center
    center_penalty_cap: float = 25.0     # max penalty

    # Bias scoring target (signed asymmetry): want positive values (antisense less stable 5' than sense)
    # If your asymmetry is signed opposite, flip with --flip-asymmetry
    bias_opt_range: Tuple[float, float] = (0.5, 4.0)
    bias_ok_range: Tuple[float, float] = (-1.0, 6.0)

    # Guide and passenger scoring weights
    lambda_guide: float = 0.50
    lambda_bias: float = 0.30
    lambda_pass: float = 0.20

    # Guide score components weights
    w_tm: float = 0.25
    w_dg: float = 0.25
    w_seed: float = 0.35
    w_end: float = 0.15

    # Passenger score components weights
    w_p_end: float = 0.45
    w_p_mod: float = 0.55

    # Seed scoring thresholds
    seed_tm_opt: float = 21.5
    seed_tm_ok: float = 24.0

    # Duplex energy "good" band (kcal/mol)
    dg_good: Tuple[float, float] = (-35.0, -20.0)
    dg_ok: Tuple[float, float] = (-40.0, -15.0)

    # Duplex Tm "good" band (°C) – подстройте под ваши данные
    tm_good: Tuple[float, float] = (48.0, 60.0)
    tm_ok: Tuple[float, float] = (45.0, 65.0)


def load_params(path: Optional[str]) -> ModParams:
    p = ModParams()
    p.zone_weights = _default_zone_weights()
    if not path:
        return p
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    # Shallow override
    for k, v in d.items():
        if hasattr(p, k):
            setattr(p, k, v)
    # Zone weights can be nested
    if "zone_weights" in d and isinstance(d["zone_weights"], dict):
        p.zone_weights = d["zone_weights"]
    return p


def recommend_baseline_pattern(sense: str, antisense: str) -> Tuple[Dict[str, Dict[str, List[int]]], Dict[str, Any]]:
    """
    Возвращает базовый паттерн модификаций (1-based позиции).
    Структура:
      pattern = {
        'guide': {'F':[...], 'OMe':[...]}   # antisense
        'passenger': {'F':[...], 'OMe':[...]}  # sense
      }

    Базовая логика (консервативная):
    - Guide (antisense):
        * OMe в seed (pos 2) как минимальный off-target "safety knob"
        * F умеренно в 3' части (после 12), через один нуклеотид
        * избегать центра 9-12 (особенно 10-11)
    - Passenger (sense):
        * более плотные OMe по чётным позициям 2-18
        * немного F ближе к 3' концу (15-19) для стабильности/PK
    """
    L = min(len(sense), len(antisense))
    # Guide
    guide_ome = [2] if L >= 2 else []
    guide_f = [i for i in range(13, L + 1, 2)]  # 13,15,17,...
    # Passenger
    pass_ome = [i for i in range(2, min(18, L) + 1, 2)]  # 2,4,6,...
    pass_f = [i for i in range(max(15, 1), L + 1, 2)]    # 15,17,19,...

    pattern = {
        "guide": {"F": guide_f, "OMe": guide_ome},
        "passenger": {"F": pass_f, "OMe": pass_ome},
    }
    meta = {"pattern_name": "baseline_v1", "L": L}
    return pattern, meta


def apply_mod_corrections(
    tm_base: float,
    dg_base: float,
    sense: str,
    antisense: str,
    pattern: Dict[str, Dict[str, List[int]]],
    params: ModParams,
) -> Dict[str, float]:
    """
    Применяет позиционно-взвешенные поправки к Tm и ΔG.
    Возвращает словарь с tm_corr_total, dg_corr_total и (опционально) раздельными коррекциями.
    """
    L = min(len(sense), len(antisense))

    def corr_for_strand(strand_name: str) -> Tuple[float, float, int]:
        # strand_name: 'guide' or 'passenger'
        zweights = params.zone_weights.get(strand_name, {})
        dT = 0.0
        dG = 0.0
        center_mods = 0
        for mod_type, positions in pattern.get(strand_name, {}).items():
            for pos1 in positions:
                if pos1 < 1 or pos1 > L:
                    continue
                zone = _zone_of_pos(pos1, L)
                w = float(zweights.get(zone, zweights.get("other", 1.0)))
                if mod_type == "F":
                    dT += w * params.alpha_F
                    dG += w * params.beta_F
                elif mod_type == "OMe":
                    dT += w * params.alpha_OMe
                    dG += w * params.beta_OMe
                # central penalty tracking only for guide
                if strand_name == "guide" and zone == "center":
                    center_mods += 1
        return dT, dG, center_mods

    dT_g, dG_g, center_mods = corr_for_strand("guide")
    dT_p, dG_p, _ = corr_for_strand("passenger")

    out = {}
    # total corr
    out["tm_corr_total"] = round(tm_base + dT_g + dT_p, 2)
    out["dg_corr_total"] = round(dg_base + dG_g + dG_p, 3)

    if params.separate_corr:
        out["tm_corr_guide"] = round(tm_base + dT_g, 2)
        out["dg_corr_guide"] = round(dg_base + dG_g, 3)
        out["tm_corr_passenger"] = round(tm_base + dT_p, 2)
        out["dg_corr_passenger"] = round(dg_base + dG_p, 3)

    # center penalty term for guide
    pen = min(params.center_penalty_cap, center_mods * params.center_penalty_per_mod)
    out["center_mods_guide"] = int(center_mods)
    out["center_penalty"] = float(pen)
    return out


# ----------------------------
# Scoring functions
# ----------------------------

def _band_score(x: float, good: Tuple[float, float], ok: Tuple[float, float]) -> float:
    """
    Преобразование значения x в 0..100 по полосам:
    - в good -> 100
    - между ok и good -> линейно
    - вне ok -> 0..40 (плавный штраф)
    """
    lo_g, hi_g = good
    lo_o, hi_o = ok
    if lo_g <= x <= hi_g:
        return 100.0
    if lo_o <= x < lo_g:
        return 100.0 * (x - lo_o) / max(1e-9, (lo_g - lo_o))
    if hi_g < x <= hi_o:
        return 100.0 * (hi_o - x) / max(1e-9, (hi_o - hi_g))
    # outside ok: soft floor
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
    """
    Предпочтение концов для strand selection:
    - guide 5' желательно A/U
    - passenger 5' желательно G/C
    """
    s = 0.0
    if guide5 in {"A", "U"}:
        s += 50.0
    if pass5 in {"G", "C"}:
        s += 50.0
    return s


def score_bias(asym_signed: Optional[float], params: ModParams) -> float:
    """
    Bias score (0..100) по signed asymmetry:
    - optimal: bias_opt_range -> 100
    - ok: bias_ok_range -> 60..90
    - outside -> down to 0..40
    """
    if asym_signed is None or (isinstance(asym_signed, float) and math.isnan(asym_signed)):
        return 50.0
    a = float(asym_signed)
    opt_lo, opt_hi = params.bias_opt_range
    ok_lo, ok_hi = params.bias_ok_range
    if opt_lo <= a <= opt_hi:
        return 100.0
    if ok_lo <= a < opt_lo:
        return 60.0 + 40.0 * (a - ok_lo) / max(1e-9, (opt_lo - ok_lo))
    if opt_hi < a <= ok_hi:
        return 60.0 + 40.0 * (ok_hi - a) / max(1e-9, (ok_hi - opt_hi))
    dist = min(abs(a - ok_lo), abs(a - ok_hi))
    return max(0.0, 40.0 - 10.0 * dist)


def score_guide(tm_corr: float, dg_corr: float, seed_tm: Optional[float],
               guide5: str, pass5: str, center_penalty: float, params: ModParams) -> float:
    tm_s = _band_score(tm_corr, params.tm_good, params.tm_ok)
    dg_s = _band_score(dg_corr, params.dg_good, params.dg_ok)
    seed_s = _seed_score(seed_tm, params.seed_tm_opt, params.seed_tm_ok)
    end_s = _end_pref_score(guide5, pass5)

    raw = (params.w_tm * tm_s + params.w_dg * dg_s + params.w_seed * seed_s + params.w_end * end_s)
    raw = max(0.0, raw - center_penalty)
    return round(min(100.0, raw), 1)


def score_passenger(pass5: str, n_mods_pass: int, L: int, params: ModParams) -> float:
    """
    Passenger score: хотим подавление passenger (т.е. чтобы passenger НЕ выбиралась guide).
    Прокси:
    - 5' passenger G/C хорошо
    - больше модификаций passenger (в разумных пределах) хорошо
    """
    end = 100.0 if pass5 in {"G", "C"} else 60.0
    # плотность модификаций
    dens = 0.0 if L <= 0 else n_mods_pass / L
    # saturation curve
    mod_s = 100.0 * (1.0 - math.exp(-4.0 * dens))  # быстро насыщается
    raw = params.w_p_end * end + params.w_p_mod * mod_s
    return round(min(100.0, raw), 1)


def count_mods(pattern: Dict[str, Dict[str, List[int]]]) -> Dict[str, int]:
    def _c(x: Dict[str, List[int]]) -> int:
        return sum(len(v) for v in x.values())
    return {
        "n_mods_guide": _c(pattern.get("guide", {})),
        "n_mods_passenger": _c(pattern.get("passenger", {})),
    }


def score_total(sg: float, sb: float, sp: float, params: ModParams) -> float:
    return round(min(100.0, params.lambda_guide * sg + params.lambda_bias * sb + params.lambda_pass * sp), 1)


# ----------------------------
# Main runner
# ----------------------------

def main():
    ap = argparse.ArgumentParser(description="Addon scoring: 2'-F vs 2'-OMe + guide/passenger scoring")
    ap.add_argument("--in", dest="inp", required=True, help="Входной CSV из ваших анализаторов")
    ap.add_argument("--out", dest="out", required=True, help="Выходной CSV (расширенный)")
    ap.add_argument("--params-json", default=None, help="JSON с параметрами модели (опционально)")
    ap.add_argument("--flip-asymmetry", action="store_true",
                    help="Инвертировать знак thermodynamic_asymmetry перед bias scoring (если ваша метрика обратная)")
    ap.add_argument("--pattern", choices=["baseline_v1"], default="baseline_v1",
                    help="Какой паттерн модификаций применять (пока реализован baseline_v1)")
    args = ap.parse_args()

    params = load_params(args.params_json)

    df = pd.read_csv(args.inp)

    required = ["sense_sequence", "antisense_sequence", "tm_celsius", "duplex_energy_kcal_mol"]
    for c in required:
        if c not in df.columns:
            raise SystemExit(f"Не найдена обязательная колонка '{c}' во входном CSV")

    # Optional columns
    if "thermodynamic_asymmetry" not in df.columns:
        df["thermodynamic_asymmetry"] = np.nan
    if "seed_tm" not in df.columns:
        df["seed_tm"] = np.nan

    rows = []
    for _, r in df.iterrows():
        sense = str(r["sense_sequence"]).upper()
        anti = str(r["antisense_sequence"]).upper()
        L = min(len(sense), len(anti))
        if L == 0:
            continue

        tm0 = float(r["tm_celsius"])
        dg0 = float(r["duplex_energy_kcal_mol"])
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

        # pattern
        pattern, meta = recommend_baseline_pattern(sense, anti)

        corr = apply_mod_corrections(tm0, dg0, sense, anti, pattern, params)
        mods = count_mods(pattern)

        # Choose which corrected values to score on:
        tm_for_guide = corr.get("tm_corr_guide", corr["tm_corr_total"])
        dg_for_guide = corr.get("dg_corr_guide", corr["dg_corr_total"])

        sg = score_guide(
            tm_corr=tm_for_guide,
            dg_corr=dg_for_guide,
            seed_tm=seed_tm_val,
            guide5=anti[0],
            pass5=sense[0],
            center_penalty=corr["center_penalty"],
            params=params
        )
        sb = score_bias(asym, params)
        sp = score_passenger(sense[0], mods["n_mods_passenger"], L, params)
        st = score_total(sg, sb, sp, params)

        out = dict(r)
        out.update(corr)
        out.update(mods)
        out["score_guide"] = sg
        out["score_bias"] = sb
        out["score_passenger"] = sp
        out["score_total"] = st
        out["mod_pattern_guide"] = json.dumps(pattern["guide"], ensure_ascii=False)
        out["mod_pattern_passenger"] = json.dumps(pattern["passenger"], ensure_ascii=False)
        out["mod_pattern_meta"] = json.dumps(meta, ensure_ascii=False)
        rows.append(out)

    outdf = pd.DataFrame(rows)

    # Sort by new score if present
    if "score_total" in outdf.columns:
        outdf = outdf.sort_values("score_total", ascending=False).reset_index(drop=True)

    outdf.to_csv(args.out, index=False, encoding="utf-8")
    print(f"✓ Saved: {args.out} ({len(outdf)} rows)")


if __name__ == "__main__":
    main()
