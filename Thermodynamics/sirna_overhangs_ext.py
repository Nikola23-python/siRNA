#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sirna_overhangs_ext.py

Расширение для siRNA_Thermo_Analyzer.py и siRNA_Thermo_Analyzer(2).py:
- добавляет поддержку оверхэнгов: dTdT, UU, и "умные" (модифицированные) оверхэнги (2'-O-Me / 2'-F).
- обеспечивает совместимость с ViennaRNA (RNA.duplexfold) даже при концептуальном dT (T -> U для энергетики).
- даёт единый интерфейс: build_overhangs(), energy_duplex_with_overhangs() и apply_overhang_corrections().

Как использовать (минимальный патч):
1) Положите этот файл рядом с вашим анализатором.
2) В siRNA_Thermo_Analyzer*.py добавьте:
     from sirna_overhangs_ext import OverhangConfig, energy_duplex_with_overhangs, apply_overhang_corrections
3) В __init__ анализатора задайте:
     self.OVERHANG = OverhangConfig(kind="dTdT")   # или "UU", "OMe_UU", "F_UU"
4) В calculate_thermodynamics замените блок расчёта energy на вызов energy_duplex_with_overhangs(...)
   и примените apply_overhang_corrections(...) при необходимости.

Примечание по "умным" оверхэнгам:
- ViennaRNA не знает о химических модификациях 2'-O-Me/2'-F. Поэтому:
  * геометрию/парность мы считаем на обычных 'U' (то есть как UU-overhang),
  * а эффект модификаций учитываем как небольшую поправку к ΔG (и опционально к Tm),
    которая включается ТОЛЬКО для оверхэнгов.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple, Dict, Any

import RNA


@dataclass
class OverhangConfig:
    """
    kind:
      - "dTdT"   : концептуально dT dT (для энергетики считаем как UU; в строке можно хранить TT)
      - "UU"     : рибонуклеотидный UU
      - "OMe_UU" : 2'-O-Me модифицированный UU (энергетику считаем как UU + корректировка)
      - "F_UU"   : 2'-F модифицированный UU (энергетику считаем как UU + корректировка)
      - "custom" : задать overhang_sense/overhang_antisense вручную (RNA-алфавит A/U/C/G)
    """
    kind: str = "dTdT"
    overhang_sense: str = "UU"
    overhang_antisense: str = "UU"

    # Empirical corrections for "smart" overhangs (very small compared to body modifications)
    # Units: (ΔTm per overhang nucleotide, ΔΔG per overhang nucleotide)
    dTm_per_nt: float = 0.0
    dG_per_nt: float = 0.0

    def resolved(self) -> "OverhangConfig":
        k = (self.kind or "").strip()
        if k == "dTdT":
            return OverhangConfig(kind="dTdT", overhang_sense="UU", overhang_antisense="UU",
                                  dTm_per_nt=0.0, dG_per_nt=0.0)
        if k == "UU":
            return OverhangConfig(kind="UU", overhang_sense="UU", overhang_antisense="UU",
                                  dTm_per_nt=0.0, dG_per_nt=0.0)
        if k == "OMe_UU":
            return OverhangConfig(kind="OMe_UU", overhang_sense="UU", overhang_antisense="UU",
                                  dTm_per_nt=0.10, dG_per_nt=-0.05)
        if k == "F_UU":
            return OverhangConfig(kind="F_UU", overhang_sense="UU", overhang_antisense="UU",
                                  dTm_per_nt=0.15, dG_per_nt=-0.08)
        if k == "custom":
            return self
        return OverhangConfig(kind="dTdT").resolved()


def build_overhangs(seq_sense: str, seq_antisense: str, cfg: OverhangConfig) -> Tuple[str, str]:
    """Добавляет 3'-оверхэнги к обеим цепям (канонически)."""
    c = cfg.resolved()
    s = (seq_sense or "").upper()
    a = (seq_antisense or "").upper()
    return s + c.overhang_sense, a + c.overhang_antisense


def build_overhangs_legacy(seq_sense: str, seq_antisense: str, cfg: OverhangConfig) -> Tuple[str, str]:
    """
    Совместимость с вашим текущим кодом:
      sense: overhang + seq
      antisense: seq + overhang
    """
    c = cfg.resolved()
    s = (seq_sense or "").upper()
    a = (seq_antisense or "").upper()
    return c.overhang_sense + s, a + c.overhang_antisense


def energy_duplex_with_overhangs(seq_sense: str, seq_antisense: str, cfg: OverhangConfig,
                                legacy_orientation: bool = True) -> float:
    """
    ΔG дуплекса (kcal/mol) с учётом оверхэнгов.
    Для OMe_UU / F_UU добавляет небольшую ΔΔG-поправку на число нуклеотидов оверхэнга.
    """
    c = cfg.resolved()
    if legacy_orientation:
        s2, a2 = build_overhangs_legacy(seq_sense, seq_antisense, c)
    else:
        s2, a2 = build_overhangs(seq_sense, seq_antisense, c)

    try:
        e = float(RNA.duplexfold(s2, a2).energy)
    except Exception:
        e = -25.0

    n_overhang = len(c.overhang_sense) + len(c.overhang_antisense)
    e_corr = e + n_overhang * float(c.dG_per_nt)
    return round(float(e_corr), 3)


def apply_overhang_corrections(tm_base: float, dg_base: float, cfg: OverhangConfig) -> Dict[str, Any]:
    """Сохраняет метаданные оверхэнгов и (опционально) применяет ΔTm-поправку для оверхэнгов."""
    c = cfg.resolved()
    n_overhang = len(c.overhang_sense) + len(c.overhang_antisense)
    tm_corr = tm_base + n_overhang * float(c.dTm_per_nt)
    return {
        "overhang_kind": c.kind,
        "overhang_sense": c.overhang_sense,
        "overhang_antisense": c.overhang_antisense,
        "overhang_nt_total": n_overhang,
        "tm_with_overhang_tm_corr": round(float(tm_corr), 2),
        "dg_with_overhang_corr": round(float(dg_base), 3),
    }
