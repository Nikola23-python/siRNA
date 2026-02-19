import json
import sys
import os
import math
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
import pandas as pd
import numpy as np

try:
    import RNA

    HAS_VIENNA = True
except ImportError:
    HAS_VIENNA = False
    print("⚠ ViennaRNA не установлен. Некоторые расчеты будут ограничены.")

try:
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq

    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    from sirna_overhangs_ext import OverhangConfig, energy_duplex_with_overhangs, apply_overhang_corrections

    HAS_OVERHANGS = True
except ImportError:
    HAS_OVERHANGS = False
    print("⚠ Модуль sirna_overhangs_ext не найден")

# Параметры модификаций 2'OMe - обновленные значения
TM_CORRECTION_PER_2OME = +0.6  # °C на каждую 2′OMe (более консервативно)
ENERGY_CORRECTION_PER_2OME = -0.25  # ккал/моль (стабилизация)


def _gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq) * 100


def _has_long_gc_stretch(seq: str, max_length: int = 7) -> bool:
    """Проверяет наличие длинных GC-секций (>7 нт) - более строго"""
    seq = seq.upper()
    pattern = r'[GC]{' + str(max_length + 1) + r',}'
    return bool(re.search(pattern, seq))


def _has_homopolymer(seq: str, max_length: int = 6) -> bool:  # Увеличить с 4 до 6
    """Проверяет наличие гомополимерных участков - менее строго"""
    seq = seq.upper()
    for base in 'AUCG':
        if base * (max_length + 1) in seq:
            return True
    return False


def _calculate_seed_tm(seq: str) -> Optional[float]:
    """Рассчитывает Tm для seed-региона (позиции 2-8, 1-based)"""
    if len(seq) < 8:
        return None
    seed_seq = seq[1:8]  # позиции 2-8 (0-based индексы 1-7)
    try:
        gc_count = seed_seq.count('G') + seed_seq.count('C')
        au_count = len(seed_seq) - gc_count
        tm_seed = 2.0 * au_count + 4.0 * gc_count  # Более точная формула
        return round(tm_seed, 1)
    except:
        return None


def _check_rna_sequence(seq: str) -> bool:
    """Проверяет корректность RNA последовательности"""
    if not seq:
        return False
    seq = seq.upper()
    valid_chars = set('AUCG')
    return all(c in valid_chars for c in seq)


def recommend_2ome_positions(sense_seq: str, antisense_seq: str):
    """
    УЛУЧШЕННЫЕ рекомендации для 2′OMe модификаций:
    - Sense: позиции 2, 4, 6, 14, 16 (1-based) для стабильности и защиты от нуклеаз
    - Antisense: только позиции 2, 14, 16 - минимум модификаций для guide
    Избегаем seed регион (2-8) на guide и 5'-концы
    """
    sense_positions = []
    antisense_positions = []

    # Sense strand (passenger) - умеренное количество модификаций
    # Позиции 2, 4, 6 - для стабильности вне seed
    for pos in [2, 4, 6]:
        if pos <= len(sense_seq):
            sense_positions.append(pos - 1)  # Convert to 0-based

    # Позиции 14, 16 - для защиты от нуклеаз в 3'-регионе
    for pos in [14, 16]:
        if pos <= len(sense_seq):
            sense_positions.append(pos - 1)

    # Antisense strand (guide) - минимальные модификации
    # Только позиции 2, 14, 16 - избегаем seed регион кроме позиции 2
    for pos in [2, 14, 16]:
        if pos <= len(antisense_seq):
            antisense_positions.append(pos - 1)

    return sense_positions, antisense_positions


class ConfigManager:
    """Управление конфигурацией анализа с РЕАЛИСТИЧНЫМИ настройками"""

    DEFAULT_CONFIG = {
        "version": "3.0",
        "overhang": {
            "kind": "dTdT",
            "options": ["dTdT", "UU", "OMe_UU", "F_UU"],
            "legacy_orientation": True
        },
        "modifications": {
            "enabled": True,
            "2OMe": True,
            "2F": False,
            "pattern_strategy": "conservative"
        },
        "scoring": {
            "weights": {
                "tm": 0.40,
                "gc": 0.15,
                "asymmetry": 0.15,
                "energy": 0.15,
                "length": 0.05,
                "seed_tm": 0.10
            },
            "thresholds": {
                "tm_optimal": 60.0,
                "tm_sigma": 5.0,
                "gc_optimal": 47.5,
                "gc_sigma": 7.5
            }
        },
        # ОСЛАБЛЕННЫЕ ФИЛЬТРЫ ДЛЯ НАЧАЛА:
        "filters": {
            "length_range": [15, 30],  # Широкий диапазон для тестирования
            "gc_range": [30, 70],      # Широкий диапазон GC
            "max_gc_stretch": 9,       # Менее строго
            "max_seed_tm": 25.0,       # Немного выше
            "min_duplex_energy": -45,  # Более широкий диапазон
            "max_duplex_energy": -15,  # Более широкий диапазон
            "asymmetry_range": [-3.0, 5.0],  # Шире
            "min_tm": 50.0,           # Ниже для начала
            "max_tm": 75.0            # Выше для начала
        },
        "quality_checks": {
            "require_correct_ends": False,  # ВЫКЛЮЧИТЬ на начальном этапе!
            "check_seed_gc": False,         # ВЫКЛЮЧИТЬ для тестирования
            "check_homopolymers": False,    # ВЫКЛЮЧИТЬ для тестирования
            "check_a_at_position_3": False, # ВЫКЛЮЧИТЬ
            "check_a_at_position_19": False # ВЫКЛЮЧИТЬ
        },
        "verbose": True  # Добавить для отладки
    }
    def __init__(self, config_file: str = None):
        self.config = self.DEFAULT_CONFIG.copy()
        if config_file and Path(config_file).exists():
            self.load(config_file)

    def load(self, filepath: str):
        with open(filepath, 'r', encoding='utf-8') as f:
            loaded = json.load(f)
            self._deep_update(self.config, loaded)

    def save(self, filepath: str):
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(self.config, f, indent=2, ensure_ascii=False)

    def _deep_update(self, target: Dict, source: Dict):
        for key, value in source.items():
            if key in target and isinstance(target[key], dict) and isinstance(value, dict):
                self._deep_update(target[key], value)
            else:
                target[key] = value

    def get_overhang_config(self) -> OverhangConfig:
        return OverhangConfig(kind=self.config["overhang"]["kind"])

    def update_setting(self, section: str, key: str, value: Any):
        keys = key.split('.')
        current = self.config[section]
        for k in keys[:-1]:
            current = current[k]
        current[keys[-1]] = value


class SiRNACoreAnalyzer:
    """Улучшенное ядро анализа siRNA"""

    def __init__(self, config: ConfigManager):
        self.config = config
        # Более реалистичные параметры для siRNA
        self.TM_PARAMS = {
            'Na': 100.0,  # Физиологическая концентрация
            'Mg': 1.0,
            'oligo_concentration': 0.2e-6,  # Типичная для siRNA
            'dNTPs': 0.0
        }

    def calculate_exact_tm_biopython(self, s: str, a: str) -> Optional[float]:
        """Точный расчет Tm с использованием Nearest-Neighbor"""
        if not BIOPYTHON_AVAILABLE:
            return None

        try:
            # Конвертируем RNA в DNA для BioPython
            s_dna = Seq(s.replace('U', 'T'))
            a_dna = Seq(a.replace('U', 'T'))

            # Расчет Tm с параметрами для siRNA
            tm_nn = mt.Tm_NN(
                s_dna,
                c_seq=a_dna,
                Na=self.TM_PARAMS['Na'],
                Mg=self.TM_PARAMS['Mg'],
                dNTPs=self.TM_PARAMS['dNTPs'],
                oligo_concentration=self.TM_PARAMS['oligo_concentration'],
                saltcorr=7  # Более подходящая коррекция для RNA
            )

            # Коррекция для RNA duplex (менее стабилен чем DNA)
            # RNA:DNA hybrid примерно на 10°C менее стабилен
            # Но siRNA duplex стабильнее из-за длины
            l = len(s)
            if l <= 21:
                corr = -4.0  # Умеренная коррекция
            elif l <= 25:
                corr = -5.0
            else:
                corr = -6.0

            # Дополнительная коррекция по GC
            gc = _gc_content(s)
            if gc < 35:
                corr -= 1.5
            elif gc > 60:
                corr += 1.5

            tm = tm_nn + corr
            return round(min(max(tm, 30), 85), 1)
        except Exception as e:
            if self.config.config.get('verbose', False):
                print(f"Ошибка расчета Tm: {e}")
            return None

    def calculate_backup_tm(self, s: str) -> float:
        """Резервный расчет Tm для siRNA (более точная формула)"""
        l = len(s)
        gc = _gc_content(s)

        # Улучшенная формула с учетом концентрации соли (100 mM Na+)
        if l <= 15:
            tm = 64.9 + 0.41 * gc - 675.0 / l
        elif l <= 20:
            tm = 69.3 + 0.41 * gc - 650.0 / l
        elif l <= 25:
            tm = 71.7 + 0.41 * gc - 625.0 / l
        else:
            tm = 73.3 + 0.41 * gc - 600.0 / l

        # Дополнительные коррекции
        if gc < 30:
            tm -= 2.5
        elif gc > 60:
            tm += 2.5

        # Коррекция для очень коротких/длинных олиго
        if l < 18:
            tm -= 1.5
        elif l > 25:
            tm += 1.5

        return round(min(max(tm, 35), 80), 1)

    def calculate_asymmetry(self, s: str, a: str) -> float:
        """Рассчитывает термодинамическую асимметрию (более точная версия)"""
        # Параметры для 4 терминальных пар
        bp_scores = {
            'GC': -2.17, 'CG': -2.17, 'AU': -1.10, 'UA': -1.10,
            'GU': -0.80, 'UG': -0.80, 'GG': -1.84, 'CC': -1.84,
            'AA': -0.90, 'UU': -0.90, 'AG': -1.28, 'GA': -1.28,
            'AC': -1.45, 'CA': -1.45, 'UC': -1.30, 'CU': -1.30
        }

        term = min(4, len(s))
        sense_stab = antisense_stab = 0.0

        # 5'-конец смысловой цепи (passenger)
        for i in range(term):
            pair = s[i] + a[-(i + 1)]
            score = bp_scores.get(pair, bp_scores.get(pair[::-1], -1.0))
            sense_stab += score

        # 5'-конец антисмысловой цепи (guide)
        for i in range(term):
            pair = a[i] + s[-(i + 1)]
            score = bp_scores.get(pair, bp_scores.get(pair[::-1], -1.0))
            antisense_stab += score

        # Асимметрия = стабильность(guide 5') - стабильность(passenger 5')
        # Положительное значение означает, что guide 5' менее стабилен
        return round(antisense_stab - sense_stab, 2)

    def calculate_thermo_score(self, tm: float, gc: float, asym: float,
                               energy: float, length: int, seed_tm: float, sense_seq: str, antisense_seq: str) -> Tuple[float, str]:
        """РЕАЛИСТИЧНЫЙ расчет комплексной оценки siRNA"""

        # 1. Tm - САМЫЙ ВАЖНЫЙ параметр (40% веса)
        if tm < 50:
            tm_s = 0  # Слишком низкий - недопустимо
        elif 58 <= tm <= 62:
            tm_s = 100  # Оптимальный диапазон
        elif 55 <= tm < 58:
            tm_s = 80 + (tm - 55) * 6.67  # 55°C=80, 58°C=100
        elif 62 < tm <= 65:
            tm_s = 100 - (tm - 62) * 6.67  # 62°C=100, 65°C=80
        elif 50 <= tm < 55:
            tm_s = 40 + (tm - 50) * 8  # 50°C=40, 55°C=80
        elif 65 < tm <= 70:
            tm_s = 80 - (tm - 65) * 4  # 65°C=80, 70°C=60
        else:
            tm_s = 20

        # 2. GC состав (15%)
        if 45 <= gc <= 50:
            gc_s = 100  # Оптимально
        elif 40 <= gc < 45:
            gc_s = 80 + (gc - 40) * 4  # 40%=80, 45%=100
        elif 50 < gc <= 55:
            gc_s = 100 - (gc - 50) * 4  # 50%=100, 55%=80
        elif 35 <= gc < 40:
            gc_s = 60 + (gc - 35) * 4  # 35%=60, 40%=80
        elif 55 < gc <= 60:
            gc_s = 80 - (gc - 55) * 4  # 55%=80, 60%=60
        else:
            gc_s = 30

        # 3. Энергия дуплекса (15%)
        if -30 <= energy <= -25:
            en_s = 100  # Идеальный диапазон
        elif -35 <= energy < -30:
            en_s = 80 + (energy + 35) * 4  # -35=80, -30=100
        elif -25 < energy <= -20:
            en_s = 100 - (energy + 20) * 4  # -20=80, -25=100
        elif -40 <= energy < -35:
            en_s = 60 + (energy + 40) * 4  # -40=60, -35=80
        elif -20 < energy <= -15:
            en_s = 80 - (energy + 15) * 4  # -15=60, -20=80
        else:
            en_s = 40

        # 4. Seed Tm (10%)
        if seed_tm is None:
            seed_s = 50
        elif 18 <= seed_tm <= 21:
            seed_s = 100  # Оптимальный диапазон для seed
        elif 15 <= seed_tm < 18:
            seed_s = 60 + (seed_tm - 15) * 13.33  # 15=60, 18=100
        elif 21 < seed_tm <= 24:
            seed_s = 100 - (seed_tm - 21) * 13.33  # 21=100, 24=60
        elif 12 <= seed_tm < 15:
            seed_s = 30 + (seed_tm - 12) * 10  # 12=30, 15=60
        elif 24 < seed_tm <= 27:
            seed_s = 60 - (seed_tm - 24) * 10  # 24=60, 27=30
        else:
            seed_s = 20

        # 5. Асимметрия (15%)
        # Желательно положительное значение (guide 5' менее стабилен)
        if 0.5 <= asym <= 3.0:
            asym_s = 100  # Оптимальный bias
        elif -0.5 <= asym < 0.5:
            asym_s = 80
        elif 3.0 < asym <= 5.0:
            asym_s = 80 - (asym - 3.0) * 10  # 3.0=80, 5.0=60
        elif -2.0 <= asym < -0.5:
            asym_s = 60 + (asym + 2.0) * 13.33  # -2.0=60, -0.5=80
        elif 5.0 < asym <= 7.0:
            asym_s = 60 - (asym - 5.0) * 10  # 5.0=60, 7.0=40
        elif -4.0 <= asym < -2.0:
            asym_s = 40 + (asym + 4.0) * 10  # -4.0=40, -2.0=60
        else:
            asym_s = 20

        # 6. Длина (5%)
        if length == 21:
            len_s = 100  # Золотой стандарт
        elif 19 <= length <= 23:
            len_s = 80
        elif 17 <= length <= 25:
            len_s = 60
        else:
            len_s = 30

        # Применяем веса
        w = self.config.config['scoring']['weights']
        score = (tm_s * w['tm'] + gc_s * w['gc'] + asym_s * w['asymmetry'] +
                 en_s * w['energy'] + len_s * w['length'] + seed_s * w['seed_tm'])

        score = round(max(0, min(100, score)), 1)

        # Категория
        if score >= 75:
            category = 'Excellent'
        elif score >= 60:
            category = 'Good'
        elif score >= 45:
            category = 'Moderate'
        else:
            category = 'Poor'

        # ДОБАВЛЯЕМ ШТРАФ ЗА НЕПРАВИЛЬНЫЕ 5'-КОНЦЫ:
        ends_penalty = 0

        # Guide (antisense) должен начинаться с A/U
        if antisense_seq[0] not in {'A', 'U'}:
            ends_penalty += 15  # Значительный штраф

            # Passenger (sense) должен начинаться с G/C
        if sense_seq[0] not in {'G', 'C'}:
            ends_penalty += 10  # Меньший штраф

            # Применяем штраф
        score = max(0, score - ends_penalty)


        return score, category

    def check_sequence_quality(self, sense_seq: str, antisense_seq: str) -> Dict[str, Any]:
        """СТРОГАЯ проверка качества siRNA пары"""
        issues = []
        warnings = []

        # Базовые проверки
        if not _check_rna_sequence(sense_seq) or not _check_rna_sequence(antisense_seq):
            return {'passed': False, 'issues': ['Некорректная последовательность'], 'warnings': []}

        length = len(sense_seq)
        if length != len(antisense_seq):
            return {'passed': False, 'issues': ['Разная длина цепей'], 'warnings': []}

        # 1. Проверка длины (критически важно)
        length_range = self.config.config['filters']['length_range']
        if not (length_range[0] <= length <= length_range[1]):
            issues.append(f"Длина {length} нт (должна быть {length_range[0]}-{length_range[1]})")
        elif length != 21:
            warnings.append(f"Длина {length} нт (оптимально 21)")

        # 2. Проверка GC состава
        gc = _gc_content(sense_seq)
        gc_range = self.config.config['filters']['gc_range']
        if not (gc_range[0] <= gc <= gc_range[1]):
            issues.append(f"GC {gc:.1f}% (должно быть {gc_range[0]}-{gc_range[1]}%)")
        elif not (45 <= gc <= 50):
            warnings.append(f"GC {gc:.1f}% (оптимально 45-50%)")

        # 3. КРИТИЧЕСКОЕ: Проверка 5'-концов
        qc_config = self.config.config['quality_checks']
        if qc_config.get('require_correct_ends', True):
            # Антисмысловая (guide) должна начинаться с A или U
            if antisense_seq[0] not in {'A', 'U'}:
                issues.append(f"5'-конец guide (antisense) должен быть A/U (у вас {antisense_seq[0]})")

            # Смысловая (passenger) должна начинаться с G или C
            if sense_seq[0] not in {'G', 'C'}:
                issues.append(f"5'-конец passenger (sense) должен быть G/C (у вас {sense_seq[0]})")
        else:
            # Только предупреждения если не требуется строго
            if antisense_seq[0] not in {'A', 'U'}:
                warnings.append(f"5'-конец guide не A/U ({antisense_seq[0]}) - может снизить эффективность")
            if sense_seq[0] not in {'G', 'C'}:
                warnings.append(f"5'-конец passenger не G/C ({sense_seq[0]}) - может снизить strand selection")

        # 4. Проверка seed региона
        if len(antisense_seq) >= 8 and qc_config.get('check_seed_gc', True):
            seed = antisense_seq[1:8]
            gc_seed = _gc_content(seed)
            if gc_seed > 65:
                issues.append(f"Seed регион слишком GC-rich: {gc_seed:.1f}% (max 65%)")
            elif gc_seed < 25:
                warnings.append(f"Seed регион AU-rich: {gc_seed:.1f}% (optimal 30-60%)")

        # 5. Проверка на гомополимеры
        if qc_config.get('check_homopolymers', True):
            if _has_homopolymer(sense_seq, 4):
                warnings.append("Гомополимеры (>4) в sense цепи")
            if _has_homopolymer(antisense_seq, 4):
                warnings.append("Гомополимеры (>4) в antisense цепи")

        # 6. Проверка длинных GC stretches
        max_gc_stretch = self.config.config['filters']['max_gc_stretch']
        if _has_long_gc_stretch(sense_seq, max_gc_stretch):
            issues.append(f"GC stretch >{max_gc_stretch} нт в sense цепи")
        if _has_long_gc_stretch(antisense_seq, max_gc_stretch):
            issues.append(f"GC stretch >{max_gc_stretch} нт в antisense цепи")

        # 7. Проверка A на позиции 3 sense цепи (Reynolds правило)
        if qc_config.get('check_a_at_position_3', True) and length >= 3:
            if sense_seq[2] != 'A':  # 1-based позиция 3
                warnings.append(f"Позиция 3 sense цепи не A ({sense_seq[2]}) - снижает эффективность")

        # 8. Проверка A на позиции 19 sense цепи
        if qc_config.get('check_a_at_position_19', True) and length >= 19:
            if sense_seq[18] != 'A':  # 1-based позиция 19
                warnings.append(f"Позиция 19 sense цепи не A ({sense_seq[18]})")

        return {
            'passed': len(issues) == 0,
            'issues': issues,
            'warnings': warnings
        }

    def analyze_pair(self, sense_seq: str, antisense_seq: str, pair_info: Dict) -> Optional[Dict]:
        """Анализ одной пары siRNA с жесткими требованиями"""
        try:
            # Быстрая проверка длины перед анализом
            if not (19 <= len(sense_seq) <= 21):
                return None

            # Проверка качества последовательности
            quality = self.check_sequence_quality(sense_seq, antisense_seq)
            if not quality['passed']:
                return None

            # Предварительный расчет Tm для фильтрации
            tm_est = self.calculate_backup_tm(sense_seq)
            tm_min = self.config.config['filters'].get('min_tm', 50.0)
            tm_max = self.config.config['filters'].get('max_tm', 70.0)

            if tm_est < tm_min or tm_est > tm_max:
                return None  # Сразу отбрасываем неподходящие по Tm

            # Точный расчет Tm
            tm_bio = self.calculate_exact_tm_biopython(sense_seq, antisense_seq)
            tm_final = tm_bio if tm_bio is not None else tm_est

            # Финальная проверка Tm
            if tm_final < 55 or tm_final > 65:
                return None  # Не в оптимальном диапазоне

            # Энергия дуплекса с оверхэнгами
            if HAS_OVERHANGS:
                overhang_cfg = self.config.get_overhang_config()
                energy = energy_duplex_with_overhangs(
                    sense_seq, antisense_seq, overhang_cfg,
                    legacy_orientation=self.config.config['overhang']['legacy_orientation']
                )
                oh_meta = apply_overhang_corrections(tm_final, energy, overhang_cfg)
            else:
                energy = -25.0
                oh_meta = {'overhang_kind': 'N/A', 'overhang_sense': '', 'overhang_antisense': ''}

            # Проверка энергии дуплекса
            min_energy = self.config.config['filters']['min_duplex_energy']
            max_energy = self.config.config['filters']['max_duplex_energy']
            if energy < min_energy or energy > max_energy:
                return None

            # Коррекции для 2'OMe
            if self.config.config['modifications']['2OMe']:
                sense_2ome, anti_2ome = recommend_2ome_positions(sense_seq, antisense_seq)
                total_2ome = len(sense_2ome) + len(anti_2ome)
                tm_corrected = round(min(tm_final + total_2ome * TM_CORRECTION_PER_2OME, 75), 1)
                energy_corrected = round(energy + total_2ome * ENERGY_CORRECTION_PER_2OME, 2)
                sense_2ome_str = str([p + 1 for p in sense_2ome])  # 1-based для пользователя
                anti_2ome_str = str([p + 1 for p in anti_2ome])  # 1-based для пользователя
            else:
                tm_corrected = round(tm_final, 1)
                energy_corrected = round(energy, 2)
                sense_2ome_str = '[]'
                anti_2ome_str = '[]'

            # GC состав
            gc_avg = (_gc_content(sense_seq) + _gc_content(antisense_seq)) / 2

            # Термодинамическая асимметрия
            asym = self.calculate_asymmetry(sense_seq, antisense_seq)

            # Seed Tm
            seed_tm = _calculate_seed_tm(antisense_seq)

            # Комплексная оценка
            score, category = self.calculate_thermo_score(
                tm_final, gc_avg, asym, energy, len(sense_seq), seed_tm
            )

            # Дополнительные проверки
            asym_range = self.config.config['filters']['asymmetry_range']
            if not (asym_range[0] <= asym <= asym_range[1]):
                quality['issues'].append(f"Асимметрия {asym:.1f} вне диапазона {asym_range}")
                quality['passed'] = False

            max_seed_tm = self.config.config['filters']['max_seed_tm']
            if seed_tm is not None and seed_tm > max_seed_tm:
                quality['warnings'].append(f"Tm seed-региона {seed_tm}°C > {max_seed_tm}°C")

            return {
                'pair_id': pair_info.get('pair_id', 'unknown'),
                'start_pos': pair_info.get('start_pos', 0),
                'sense_sequence': sense_seq,
                'antisense_sequence': antisense_seq,
                'original_length': len(sense_seq),
                'tm_celsius': round(tm_final, 1),
                'gc_content_percent': round(gc_avg, 1),
                'thermodynamic_asymmetry': asym,
                'duplex_energy_kcal_mol': round(energy, 2),
                'seed_tm': seed_tm,
                'thermo_score': score,
                'thermo_category': category,
                'quality_check_passed': quality['passed'],
                'quality_issues': '; '.join(quality['issues']),
                'quality_warnings': '; '.join(quality['warnings']),
                'sense_2ome_positions': sense_2ome_str,
                'antisense_2ome_positions': anti_2ome_str,
                'tm_corrected_for_2ome': tm_corrected,
                'energy_corrected_for_2ome': energy_corrected,
                'overhang_kind': oh_meta.get('overhang_kind', 'N/A'),
                'overhang_sense': oh_meta.get('overhang_sense', ''),
                'overhang_antisense': oh_meta.get('overhang_antisense', '')
            }
        except Exception as e:
            if self.config.config.get('verbose', False):
                print(f"Ошибка анализа пары {pair_info.get('pair_id', 'unknown')}: {e}")
            return None


class SiRNAInteractiveAnalyzer:
    """Интерактивный анализатор siRNA (упрощенный, только CSV)"""

    def __init__(self, config: ConfigManager = None):
        self.config = config or ConfigManager()
        self.data = None
        self.results = None
        self.core_analyzer = SiRNACoreAnalyzer(self.config)
        self.max_pairs_to_analyze = None  # Лимит для тестирования

    def interactive_menu(self):
        """Главное меню"""
        while True:
            print("\n" + "=" * 60)
            print("siRNA ACCURATE ANALYZER v3.0")
            print("=" * 60)
            print("1) Загрузить последовательности из CSV")
            print("2) Настроить параметры анализа")
            print("3) Запустить анализ кандидатов")
            print("4) Экспорт результатов в CSV")
            print("5) Сохранить/загрузить конфигурацию")
            print("6) Выход")
            print("=" * 60)

            choice = input("\nВыберите действие [1-6]: ").strip()

            if choice == "1":
                self.load_sequences()
            elif choice == "2":
                self.configure_analysis()
            elif choice == "3":
                self.run_analysis()
            elif choice == "4":
                self.export_results()
            elif choice == "5":
                self.config_menu()
            elif choice == "6":
                print("\nВыход из программы. До свидания!")
                break
            else:
                print("Неверный выбор. Попробуйте снова.")

    def load_sequences(self):
        """Загрузка последовательностей из CSV файла"""
        print("\n" + "-" * 60)
        print("ЗАГРУЗКА ПОСЛЕДОВАТЕЛЬНОСТЕЙ ИЗ CSV")
        print("-" * 60)

        filepath = input("Путь к CSV файлу: ").strip()
        if not filepath:
            print("❌ Не указан путь к файлу")
            return

        if not Path(filepath).exists():
            print(f"❌ Файл не найден: {filepath}")
            return

        try:
            self.data = pd.read_csv(filepath)

            # Проверяем наличие нужных колонок
            if 'sequence' in self.data.columns:
                print(f"✓ Загружено {len(self.data)} последовательностей")
                # Убираем пробелы и преобразуем в верхний регистр
                self.data['sequence'] = self.data['sequence'].astype(str).str.strip().str.upper()
                print(f"  Пример: {self.data['sequence'].iloc[0][:50]}...")

                # Проверяем качество последовательностей
                valid_count = sum(1 for seq in self.data['sequence'] if _check_rna_sequence(seq))
                print(f"  Корректных RNA последовательностей: {valid_count}/{len(self.data)}")

            elif 'sense_sequence' in self.data.columns and 'antisense_sequence' in self.data.columns:
                print(f"✓ Загружено {len(self.data)} пар siRNA")
                # Очищаем и проверяем обе колонки
                self.data['sense_sequence'] = self.data['sense_sequence'].astype(str).str.strip().str.upper()
                self.data['antisense_sequence'] = self.data['antisense_sequence'].astype(str).str.strip().str.upper()
                self.data['sequence'] = self.data['sense_sequence']  # Для совместимости
            else:
                print("❌ Файл должен содержать колонку 'sequence' или 'sense_sequence'/'antisense_sequence'")
                print("   Доступные колонки:", list(self.data.columns))
                self.data = None

        except Exception as e:
            print(f"❌ Ошибка загрузки CSV: {e}")
            self.data = None

    def configure_analysis(self):
        """Настройка параметров анализа"""
        print("\n" + "-" * 60)
        print("НАСТРОЙКА ПАРАМЕТРОВ АНАЛИЗА")
        print("-" * 60)

        while True:
            print("\nТекущие настройки:")
            print(f"1) Тип оверхэнга: {self.config.config['overhang']['kind']}")
            print(f"2) Модификации 2'OMe: {self.config.config['modifications']['2OMe']}")
            print(f"3) Диапазон длины: {self.config.config['filters']['length_range']}")
            print(f"4) Диапазон GC: {self.config.config['filters']['gc_range']}")
            print(
                f"5) Диапазон Tm: {self.config.config['filters'].get('min_tm', 55.0)}-{self.config.config['filters'].get('max_tm', 65.0)}°C")
            print(f"6) Максимальное количество пар для анализа: {self.max_pairs_to_analyze or 'все'}")
            print("7) Возврат в меню")

            choice = input("\nЧто изменить? [1-7]: ").strip()

            if choice == "1":
                print("Доступные типы: dTdT, UU, OMe_UU, F_UU")
                new_kind = input("Новый тип: ").strip()
                if new_kind in ["dTdT", "UU", "OMe_UU", "F_UU"]:
                    self.config.update_setting("overhang", "kind", new_kind)
                    print(f"✓ Тип оверхэнга изменен на {new_kind}")
                else:
                    print("❌ Неподдерживаемый тип")

            elif choice == "2":
                current = self.config.config['modifications']['2OMe']
                new_val = not current
                self.config.update_setting("modifications", "2OMe", new_val)
                print(f"✓ 2'OMe модификации: {'включены' if new_val else 'выключены'}")

            elif choice == "3":
                print(f"Текущий диапазон: {self.config.config['filters']['length_range']}")
                print("Рекомендуется: 19-21 (21 - оптимально)")
                min_len = input("Минимальная длина (по умолчанию 19): ").strip()
                max_len = input("Максимальная длина (по умолчанию 21): ").strip()
                try:
                    new_range = [int(min_len) if min_len else 19,
                                 int(max_len) if max_len else 21]
                    if new_range[0] < 15 or new_range[1] > 30:
                        print("⚠ Диапазон выходит за рекомендуемые пределы (15-30)")
                    self.config.update_setting("filters", "length_range", new_range)
                    print(f"✓ Новый диапазон: {new_range}")
                except ValueError:
                    print("❌ Некорректные значения")

            elif choice == "4":
                print(f"Текущий диапазон: {self.config.config['filters']['gc_range']}")
                print("Рекомендуется: 40-55 (45-50 оптимально)")
                min_gc = input("Минимальный GC% (по умолчанию 40): ").strip()
                max_gc = input("Максимальный GC% (по умолчанию 55): ").strip()
                try:
                    new_range = [float(min_gc) if min_gc else 40.0,
                                 float(max_gc) if max_gc else 55.0]
                    if new_range[0] < 20 or new_range[1] > 80:
                        print("⚠ Диапазон выходит за разумные пределы (20-80)")
                    self.config.update_setting("filters", "gc_range", new_range)
                    print(f"✓ Новый диапазон: {new_range}")
                except ValueError:
                    print("❌ Некорректные значения")

            elif choice == "5":
                min_tm = self.config.config['filters'].get('min_tm', 55.0)
                max_tm = self.config.config['filters'].get('max_tm', 65.0)
                print(f"Текущий диапазон: {min_tm}-{max_tm}°C")
                print("Рекомендуется: 55-65°C (58-62 оптимально)")
                new_min = input(f"Минимальный Tm (°C) [по умолчанию {min_tm}]: ").strip()
                new_max = input(f"Максимальный Tm (°C) [по умолчанию {max_tm}]: ").strip()
                try:
                    if new_min:
                        self.config.config['filters']['min_tm'] = float(new_min)
                    if new_max:
                        self.config.config['filters']['max_tm'] = float(new_max)
                    print(
                        f"✓ Новый диапазон Tm: {self.config.config['filters']['min_tm']}-{self.config.config['filters']['max_tm']}°C")
                except ValueError:
                    print("❌ Некорректные значения")

            elif choice == "6":
                print(f"Текущее значение: {self.max_pairs_to_analyze or 'нет ограничения'}")
                print("Установите лимит для тестирования или 0 для всех")
                new_val = input("Макс. количество пар: ").strip()
                try:
                    val = int(new_val)
                    self.max_pairs_to_analyze = val if val > 0 else None
                    print(f"✓ Макс. количество пар: {self.max_pairs_to_analyze or 'все'}")
                except ValueError:
                    print("❌ Некорректное значение")

            elif choice == "7":
                break


    def create_complementary_pairs(self) -> pd.DataFrame:
        """Создает комплементарные пары из sense последовательностей"""
        if self.data is None or self.data.empty:
            return pd.DataFrame()

        comp_dict = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        pairs = []

        print("\nСоздание комплементарных пар...")

        # Ограничиваем количество если задано
        data_to_process = self.data
        if self.max_pairs_to_analyze and len(self.data) > self.max_pairs_to_analyze:
            data_to_process = self.data.head(self.max_pairs_to_analyze)
            print(f"Ограничение: обработано первых {self.max_pairs_to_analyze} последовательностей")

        # Отладочные счетчики
        counters = {
            'total': 0,
            'invalid_seq': 0,
            'length_filter': 0,
            'gc_filter': 0,
            'homopolymer_filter': 0,
            'ends_filter': 0,
            'passed': 0
        }

        for idx, row in data_to_process.iterrows():
            counters['total'] += 1
            seq = str(row['sequence']).upper()

            # Проверка корректности последовательности
            if not _check_rna_sequence(seq):
                counters['invalid_seq'] += 1
                continue

            # Быстрые фильтры перед созданием пары
            length = len(seq)
            length_range = self.config.config['filters']['length_range']
            if not (length_range[0] <= length <= length_range[1]):
                counters['length_filter'] += 1
                continue

            # Проверка GC
            gc_val = _gc_content(seq)
            gc_range = self.config.config['filters']['gc_range']
            if not (gc_range[0] <= gc_val <= gc_range[1]):
                counters['gc_filter'] += 1
                continue

            # Проверка на гомополимеры
            if _has_homopolymer(seq, 5):  # Обратите внимание: 5 - слишком строго!
                counters['homopolymer_filter'] += 1
                continue

            # Создаем антисмысловую цепь
            anti = ''.join(comp_dict.get(n, n) for n in seq)

            # Быстрая проверка 5'-концов
            if (self.config.config['quality_checks']['require_correct_ends'] and
                    (anti[0] not in {'A', 'U'} or seq[0] not in {'G', 'C'})):
                counters['ends_filter'] += 1
                continue

            pairs.append({
                'pair_id': f"pair_{idx:06d}",
                'sense_sequence': seq,
                'antisense_sequence': anti,
                'original_length': length,
                'start_pos': idx,
                'gc_content_percent': round(gc_val, 1)
            })
            counters['passed'] += 1

        print(f"\nСтатистика фильтрации:")
        print(f"  Всего проверено: {counters['total']}")
        print(f"  Некорректные последовательности: {counters['invalid_seq']}")
        print(f"  Отфильтровано по длине: {counters['length_filter']}")
        print(f"  Отфильтровано по GC: {counters['gc_filter']}")
        print(f"  Отфильтровано по гомополимерам: {counters['homopolymer_filter']}")
        print(f"  Отфильтровано по 5'-концам: {counters['ends_filter']}")
        print(f"  Прошло фильтры: {counters['passed']}")

        if counters['passed'] == 0:
            print("\n⚠ Все последовательности отфильтрованы!")
            print("Основные причины:")
            if counters['gc_filter'] > counters['total'] * 0.5:
                print(
                    f"  • GC состав: {counters['gc_filter']}/{counters['total']} не в диапазоне {self.config.config['filters']['gc_range']}")
            if counters['ends_filter'] > counters['total'] * 0.5:
                print(f"  • 5'-концы: {counters['ends_filter']}/{counters['total']} не соответствуют требованиям")
                print(f"    (sense должен начинаться с G/C, antisense с A/U)")
            if counters['homopolymer_filter'] > counters['total'] * 0.5:
                print(
                    f"  • Гомополимеры: {counters['homopolymer_filter']}/{counters['total']} содержат гомополимеры >5 нт")

        print(f"✓ Создано {len(pairs)} пар для анализа")
        return pd.DataFrame(pairs)

    def run_analysis(self):
        """Запуск анализа с текущими настройками"""
        if self.data is None or self.data.empty:
            print("❌ Нет данных для анализа. Сначала загрузите последовательности.")
            return

        print("\n" + "-" * 60)
        print("ЗАПУСК АНАЛИЗА siRNA КАНДИДАТОВ")
        print("-" * 60)

        print("\nКонфигурация анализа:")
        print(f"  • Тип оверхэнга: {self.config.config['overhang']['kind']}")
        print(f"  • 2'OMe модификации: {self.config.config['modifications']['2OMe']}")
        print(f"  • Диапазон длины: {self.config.config['filters']['length_range']}")
        print(f"  • Диапазон GC: {self.config.config['filters']['gc_range']}")
        print(
            f"  • Диапазон Tm: {self.config.config['filters'].get('min_tm', 55.0)}-{self.config.config['filters'].get('max_tm', 65.0)}°C")

        # Создаем пары
        pairs_df = self.create_complementary_pairs()

        if pairs_df.empty:
            print("❌ Нет подходящих пар после первичной фильтрации")
            print("Попробуйте ослабить фильтры в настройках")
            return

        # Анализируем пары
        print(f"\nАнализ {len(pairs_df)} пар...")

        results = []
        start_time = datetime.now()

        for idx, row in pairs_df.iterrows():
            # Прогресс каждые 100 пар или 10%
            if idx > 0 and (idx % 100 == 0 or idx == len(pairs_df) - 1):
                elapsed = (datetime.now() - start_time).total_seconds()
                percent = (idx + 1) / len(pairs_df) * 100
                print(f"  Обработано {idx + 1}/{len(pairs_df)} пар ({percent:.1f}%)...")

            result = self.core_analyzer.analyze_pair(
                row['sense_sequence'],
                row['antisense_sequence'],
                row.to_dict()
            )

            if result:
                results.append(result)

        elapsed_time = (datetime.now() - start_time).total_seconds()

        if results:
            self.results = pd.DataFrame(results)
            self.results = self.results.sort_values('thermo_score', ascending=False).reset_index(drop=True)

            print(f"\n{'=' * 60}")
            print("АНАЛИЗ ЗАВЕРШЕН!")
            print('=' * 60)
            print(f"  Время анализа: {elapsed_time:.1f} сек")
            print(f"  Всего проанализировано пар: {len(pairs_df)}")
            print(f"  Прошли все фильтры: {len(self.results)}")

            if len(self.results) > 0:
                print(f"  Процент прошедших: {len(self.results) / len(pairs_df) * 100:.1f}%")

                # Статистика по категориям
                categories = self.results['thermo_category'].value_counts()
                print("\n  Распределение по категориям:")
                for cat, count in categories.items():
                    print(f"    {cat}: {count} ({count / len(self.results) * 100:.1f}%)")

                # Средние показатели
                print("\n  Средние значения:")
                print(f"    Tm: {self.results['tm_celsius'].mean():.1f} ± {self.results['tm_celsius'].std():.1f}°C")
                print(
                    f"    GC: {self.results['gc_content_percent'].mean():.1f} ± {self.results['gc_content_percent'].std():.1f}%")
                print(
                    f"    Score: {self.results['thermo_score'].mean():.1f} ± {self.results['thermo_score'].std():.1f}")

                # Топ кандидаты
                top_n = min(10, len(self.results))
                print(f"\n  Топ-{top_n} кандидатов:")
                for i in range(top_n):
                    cand = self.results.iloc[i]
                    issues = cand['quality_issues']
                    warnings = cand['quality_warnings']

                    print(f"\n  {i + 1}. {cand['pair_id']} [{cand['thermo_category']}]")
                    print(f"     Score: {cand['thermo_score']:.1f}, Tm: {cand['tm_celsius']}°C")
                    print(f"     GC: {cand['gc_content_percent']:.1f}%, Длина: {cand['original_length']} нт")
                    print(f"     Sense: {cand['sense_sequence']}")
                    print(f"     Anti:  {cand['antisense_sequence']}")

                    if warnings:
                        warnings_list = warnings.split('; ')
                        if len(warnings_list) > 2:
                            print(f"     ⚠ Предупреждения: {warnings_list[0]}; {warnings_list[1]}; ...")
                        else:
                            print(f"     ⚠ Предупреждения: {warnings}")
            else:
                print("\n⚠ Нет кандидатов, прошедших все фильтры")
                print("Попробуйте:")
                print("  • Расширить диапазоны фильтров")
                print("  • Проверить качество исходных последовательностей")
                print("  • Изменить требования к 5'-концам")
        else:
            print("❌ Анализ не дал результатов")

    def export_results(self):
        """Экспорт результатов в CSV"""
        if self.results is None or self.results.empty:
            print("❌ Нет результатов для экспорта.")
            return

        print("\n" + "-" * 60)
        print("ЭКСПОРТ РЕЗУЛЬТАТОВ В CSV")
        print("-" * 60)

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        default_name = f"sirna_analysis_{timestamp}.csv"

        filename = input(f"Имя файла [по умолчанию: {default_name}]: ").strip()
        if not filename:
            filename = default_name

        # Добавляем расширение .csv если нет
        if not filename.endswith('.csv'):
            filename += '.csv'

        output_dir = "./results"
        Path(output_dir).mkdir(exist_ok=True)

        filepath = f"{output_dir}/{filename}"

        try:
            # Сохраняем в CSV
            self.results.to_csv(filepath, index=False, encoding='utf-8')

            print(f"\n✓ Результаты сохранены: {filepath}")
            print(f"  Записей: {len(self.results)}")
            print(f"  Размер файла: {Path(filepath).stat().st_size / 1024:.1f} KB")

            # Статистика по файлу
            if len(self.results) > 0:
                print(f"\n  Категории в файле:")
                categories = self.results['thermo_category'].value_counts()
                for cat, count in categories.items():
                    print(f"    {cat}: {count}")

                # Информация о лучшем кандидате
                best = self.results.iloc[0]
                print(f"\n  Лучший кандидат в файле:")
                print(f"    ID: {best['pair_id']}, Score: {best['thermo_score']}, Tm: {best['tm_celsius']}°C")

        except Exception as e:
            print(f"❌ Ошибка при сохранении CSV: {e}")

    def config_menu(self):
        """Меню управления конфигурацией"""
        print("\n" + "-" * 60)
        print("УПРАВЛЕНИЕ КОНФИГУРАЦИЕЙ")
        print("-" * 60)

        while True:
            print("\n1) Сохранить текущую конфигурацию")
            print("2) Загрузить конфигурацию из файла")
            print("3) Сбросить к настройкам по умолчанию")
            print("4) Показать текущие настройки")
            print("5) Назад")

            choice = input("\nВыберите действие [1-5]: ").strip()

            if choice == "1":
                filename = input("Имя файла конфигурации [config.json]: ").strip()
                if not filename:
                    filename = "config.json"
                if not filename.endswith('.json'):
                    filename += '.json'
                self.config.save(filename)
                print(f"✓ Конфигурация сохранена в {filename}")

            elif choice == "2":
                filename = input("Имя файла: ").strip()
                if Path(filename).exists():
                    self.config.load(filename)
                    print(f"✓ Конфигурация загружена из {filename}")
                else:
                    print(f"❌ Файл {filename} не найден")

            elif choice == "3":
                self.config = ConfigManager()
                print("✓ Настройки сброшены к значениям по умолчанию")

            elif choice == "4":
                print("\nТекущая конфигурация:")
                print(json.dumps(self.config.config, indent=2, ensure_ascii=False))

            elif choice == "5":
                break


def main():
    """Основная функция"""
    print("\n" + "=" * 70)
    print("siRNA ACCURATE ANALYZER v3.0")
    print("=" * 70)
    print("Улучшенная система анализа siRNA с повышенной точностью")
    print("Только работа с CSV файлами")
    print("=" * 70)

    # Проверка зависимостей
    if not HAS_VIENNA:
        print("⚠ Внимание: ViennaRNA не установлен")
        print("  Расчет энергии дуплекса будет ограничен")
        print("  Установите: pip install ViennaRNA")

    if not BIOPYTHON_AVAILABLE:
        print("⚠ Внимание: BioPython не установлен")
        print("  Расчет Tm будет использовать упрощенные формулы")
        print("  Установите: pip install biopython")

    if not HAS_OVERHANGS:
        print("⚠ Внимание: модуль sirna_overhangs_ext не найден")
        print("  Функции оверхэнгов будут ограничены")
        print("  Убедитесь, что файл sirna_overhangs_ext.py находится в той же папке")

    # Создание анализатора
    analyzer = SiRNAInteractiveAnalyzer()

    # Запуск интерактивного меню
    analyzer.interactive_menu()


if __name__ == "__main__":
    main()