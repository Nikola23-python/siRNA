from datetime import datetime
import RNA
import pandas as pd
import numpy as np
import glob
import warnings
from prepare_rna import df_sense, df_antisense
from helper import helper

# Отключаем DtypeWarning
warnings.filterwarnings('ignore', category=pd.errors.DtypeWarning)


class analyz_rna:
    def __init__(self, merged_df=None):
        if merged_df is not None:
            # Если передан объединенный датафрейм, используем его
            self.df_sense = merged_df[['fragment_id', 'sense_sequence']].rename(
                columns={'sense_sequence': 'sequence'}
            )
            self.df_antisense = merged_df[['fragment_id', 'antisense_sequence']].rename(
                columns={'antisense_sequence': 'sequence'}
            )
            self.merged_df = merged_df
            print(f"📦 Загружено {len(self.df_sense)} последовательностей из объединенного датафрейма")
        else:
            # Если нет, используем исходные prepare_rna
            self.df_sense = df_sense
            self.df_antisense = df_antisense
            self.merged_df = None
            print(f"📦 Используются исходные последовательности из prepare_rna")

    """Обработка sense-последовательности"""

    def analyze_gc(self):
        """GC content of 36–52%"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            g_count = seq_upper.count('G')
            c_count = seq_upper.count('C')
            gc_count = g_count + c_count

            if length > 0:
                gc_percent = (gc_count / length) * 100
            else:
                gc_percent = 0.0

            gc_point = 0 if gc_percent < 36 else 0 if gc_percent > 52 else 1

            results.append({
                'gc_percent': round(gc_percent, 2),
                'gc_point': gc_point,
            })
        return pd.DataFrame(results)

    def u_at_10(self):
        """U at position 10 of sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 10:
                has_u_at_10 = (seq_upper[9] == 'U')
                u_tenth_point = 1 if has_u_at_10 else 0
            else:
                u_tenth_point = 0
            results.append({'u_10_point': u_tenth_point})
        return pd.DataFrame(results)

    def no_g_at_13(self):
        """No G at position 13 of sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 13:
                hasnt_g_at_13 = (seq_upper[12] != 'G')
                no_g_thirteenth_point = 1 if hasnt_g_at_13 else 0
            else:
                no_g_thirteenth_point = 0
            results.append({'no_g_13_point': no_g_thirteenth_point})
        return pd.DataFrame(results)

    def no_g_c_at_19(self):
        """No G/C at position 19 of sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 19:
                nineteenth = seq_upper[18]
                no_g_c = (nineteenth != 'G') and (nineteenth != 'C')
                no_g_c_point = 1 if no_g_c else 0
            else:
                no_g_c_point = 0
            results.append({'no_g_c_19_point': no_g_c_point})
        return pd.DataFrame(results)

    def a_at_3_19(self):
        """A at positions 3 and 19 of sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 19:
                a_at_positions = (seq_upper[2] == 'A') and (seq_upper[18] == 'A')
                a_point = 1 if a_at_positions else 0
            else:
                a_point = 0
            results.append({'a_at_3_19_point': a_point})
        return pd.DataFrame(results)

    def strong_pairing_5_end(self):
        """G/C at 5'-end of sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 1:
                g_c_5_end = (seq_upper[0] == 'G') or (seq_upper[0] == 'C')
                g_c_point = 1 if g_c_5_end else 0
            else:
                g_c_point = 0
            results.append({'g_c_5_end_point': g_c_point})
        return pd.DataFrame(results)

    def no_u_at_position_1(self):
        """Position 1 of sense strand should not be U"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 1:
                not_u = seq_upper[0] != 'U'
                pos1_point = 1 if not_u else 0
            else:
                pos1_point = 0
            results.append({'pos1_not_u_point': pos1_point})
        return pd.DataFrame(results)

    def au_at_15_19(self):
        """At least 3 A/U in positions 15-19 of sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 19:
                region = seq_upper[14:19]
                au_count = sum(1 for base in region if base in ['A', 'U'])
                au_point = 1 if au_count >= 3 else 0
            else:
                au_point = 0
            results.append({'au_15_19_point': au_point})
        return pd.DataFrame(results)

    def frequency_repeats(self):
        """Check repeats: GC ≤ 3 in a row, AT ≤ 4 in a row"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()

            max_gc_repeat = 0
            max_at_repeat = 0
            current_repeat = 1

            for i in range(1, len(seq_upper)):
                if seq_upper[i] == seq_upper[i - 1]:
                    current_repeat += 1
                else:
                    if seq_upper[i - 1] in ['G', 'C']:
                        max_gc_repeat = max(max_gc_repeat, current_repeat)
                    elif seq_upper[i - 1] in ['A', 'T', 'U']:
                        max_at_repeat = max(max_at_repeat, current_repeat)
                    current_repeat = 1

            # Проверка последнего повтора
            if seq_upper[-1] in ['G', 'C']:
                max_gc_repeat = max(max_gc_repeat, current_repeat)
            elif seq_upper[-1] in ['A', 'T', 'U']:
                max_at_repeat = max(max_at_repeat, current_repeat)

            repeats_ok = (max_gc_repeat <= 3) and (max_at_repeat <= 4)
            repeats_point = 1 if repeats_ok else 0

            results.append({'repeats_point': repeats_point})
        return pd.DataFrame(results)

    """Обработка antisense-последовательности"""

    def a_at_6(self):
        """A at position 6 of antisense strand"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 6:
                a_at_6 = (seq_upper[5] == 'A')
                a_at_6_point = 1 if a_at_6 else 0
            else:
                a_at_6_point = 0
            results.append({'a_at_6_point': a_at_6_point})
        return pd.DataFrame(results)

    def weak_pairing_5_end(self):
        """A/U at 5'-end of antisense strand"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 1:
                a_u_5_end = (seq_upper[0] == 'A') or (seq_upper[0] == 'U')
                a_u_point = 1 if a_u_5_end else 0
            else:
                a_u_point = 0
            results.append({'a_u_5_end_point': a_u_point})
        return pd.DataFrame(results)

    def gc_profile_antisense(self):
        """GC profile: 19% for positions 2-7, 52% for positions 8-18"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()

            if len(seq_upper) >= 18:
                seed_region = seq_upper[1:7]
                seed_gc = (seed_region.count('G') + seed_region.count('C')) / len(seed_region) * 100

                central_region = seq_upper[7:18]
                central_gc = (central_region.count('G') + central_region.count('C')) / len(central_region) * 100

                seed_ok = 15 <= seed_gc <= 25
                central_ok = 45 <= central_gc <= 60
                gc_profile_point = 2 if (seed_ok and central_ok) else 1 if (seed_ok or central_ok) else 0
            else:
                gc_profile_point = 0

            results.append({'gc_profile_point': gc_profile_point})
        return pd.DataFrame(results)

    def a_at_10_antisense(self):
        """A at position 10 of antisense strand"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 10:
                a_at_10 = seq_upper[9] == 'A'
                a10_point = 1 if a_at_10 else 0
            else:
                a10_point = 0
            results.append({'a10_antisense_point': a10_point})
        return pd.DataFrame(results)

    def no_a_at_19_antisense(self):
        """Position 19 of antisense strand should not be A"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()
            if len(seq_upper) >= 19:
                not_a = seq_upper[18] != 'A'
                a19_point = 1 if not_a else 0
            else:
                a19_point = 0
            results.append({'pos19_not_a_point': a19_point})
        return pd.DataFrame(results)

    """Объединение результатов"""

    def get_combined_data(self):
        """Получение объединенного датафрейма со всеми признаками"""

        # Все методы анализа
        df1 = self.analyze_gc().add_prefix('m1_')
        df2 = self.frequency_repeats().add_prefix('m2_')
        df3 = self.u_at_10().add_prefix('m3_')
        df4 = self.no_g_at_13().add_prefix('m4_')
        df5 = self.no_g_c_at_19().add_prefix('m5_')
        df6 = self.a_at_3_19().add_prefix('m6_')
        df7 = self.strong_pairing_5_end().add_prefix('m7_')
        df8 = self.no_u_at_position_1().add_prefix('m8_')
        df9 = self.au_at_15_19().add_prefix('m9_')
        df10 = self.a_at_6().add_prefix('m10_')
        df11 = self.weak_pairing_5_end().add_prefix('m11_')
        df12 = self.gc_profile_antisense().add_prefix('m12_')
        df13 = self.a_at_10_antisense().add_prefix('m13_')
        df14 = self.no_a_at_19_antisense().add_prefix('m14_')

        # Объединяем все датафреймы
        combined_df = pd.concat([
            df1, df2, df3, df4, df5, df6, df7, df8,
            df9, df10, df11, df12, df13, df14
        ], axis=1)

        # Собираем все колонки с баллами
        point_columns = [col for col in combined_df.columns if 'point' in col]
        combined_df['design_score'] = combined_df[point_columns].sum(axis=1)
        combined_df['sense_sequence'] = self.df_sense['sequence'].values
        combined_df['antisense_sequence'] = self.df_antisense['sequence'].values

        # Добавляем fragment_id из исходных данных
        combined_df['fragment_id'] = self.df_sense['fragment_id'].values

        # Если есть объединенный датафрейм, добавляем BLAST и SNP баллы
        if self.merged_df is not None:
            combined_df['blast_score'] = self.merged_df['blast_score'].values
            combined_df['snp_score'] = self.merged_df['snp_score'].values
            combined_df['total_score'] = self.merged_df['total_score'].values
            combined_df['category'] = self.merged_df['category'].values

            # Создаем финальный общий балл
            combined_df['final_score'] = (
                    combined_df['design_score'] +
                    combined_df['blast_score'] +
                    combined_df['snp_score']
            )

            # Категория на основе всех анализов
            def get_final_category(row):
                if row['final_score'] >= 20:
                    return 'Excellent'
                elif row['final_score'] >= 15:
                    return 'Good'
                elif row['final_score'] >= 10:
                    return 'Average'
                else:
                    return 'Poor'

            combined_df['final_category'] = combined_df.apply(get_final_category, axis=1)

        return combined_df


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("АНАЛИЗ EXCELLENT siRNA ИЗ MERGE_BLAST_SNP")
    print("=" * 60)

    try:
        # Ищем файл с отличными siRNA
        excellent_files = glob.glob('merged_sirna_excellent_*.csv')

        if not excellent_files:
            print("❌ Не найдены файлы merged_sirna_excellent_*.csv")
            print("   Сначала запустите merge_blast_snp.py")
            exit()

        # Берем самый свежий файл
        latest_excellent = max(excellent_files)
        print(f"📂 Загрузка: {latest_excellent}")

        # Загружаем с правильными типами данных
        dtype_dict = {
            'fragment_id': 'string',
            'sense_sequence': 'string',
            'antisense_sequence': 'string',
            'blast_score': 'Int64',
            'snp_score': 'Int64',
            'total_score': 'Int64',
            'category': 'string'
        }

        excellent_df = pd.read_csv(latest_excellent, dtype=dtype_dict, low_memory=False)
        print(f"✅ Загружено {len(excellent_df)} excellent siRNA")

        # Анализируем
        analyzer = analyz_rna(excellent_df)
        results = analyzer.get_combined_data()

        # Выводим статистику
        print(f"\n📊 СТАТИСТИКА АНАЛИЗА:")
        print(f"   Всего проанализировано: {len(results)} siRNA")
        print(f"   Средний design_score: {results['design_score'].mean():.2f}")
        print(f"   Мин design_score: {results['design_score'].min()}")
        print(f"   Макс design_score: {results['design_score'].max()}")
        print(f"   Средний final_score: {results['final_score'].mean():.2f}")

        # Распределение по финальным категориям
        print(f"\n🏷️ ФИНАЛЬНЫЕ КАТЕГОРИИ:")
        cat_stats = results['final_category'].value_counts()
        for cat, count in cat_stats.items():
            print(f"   {cat}: {count} ({count / len(results) * 100:.1f}%)")

        # Топ-10 лучших
        print(f"\n🏆 ТОП-10 ЛУЧШИХ siRNA:")
        top10 = results.nlargest(10, 'final_score')
        for idx, row in top10.iterrows():
            print(f"   {row['fragment_id']}: final_score={row['final_score']} "
                  f"(design={row['design_score']}, blast={row['blast_score']}, snp={row['snp_score']})")

        # Сохраняем результаты
        output_file = f"final_excellent_analysis_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv"
        results.to_csv(output_file, index=False)
        print(f"\n💾 Результаты сохранены в: {output_file}")

        # Сохраняем только лучшие из лучших (final_score >= 18)
        best_of_best = results[results['final_score'] >= 18]
        if len(best_of_best) > 0:
            best_file = f"best_of_excellent_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv"
            best_of_best.to_csv(best_file, index=False)
            print(f"💾 Лучшие из лучших (score≥18): {best_file} ({len(best_of_best)} siRNA)")

    except Exception as e:
        print(f"❌ Ошибка: {e}")
        import traceback

        traceback.print_exc()