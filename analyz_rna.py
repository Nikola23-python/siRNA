from datetime import datetime

import pandas as pd

from prepare_rna import df_sense, df_antisense
from helper import helper

class analyz_rna:
    def __init__(self):
        self.df_sense = df_sense
        self.df_antisense = df_antisense

    """Обработка sense-последовательности"""

    def analyze_gc(self):
        """GC content of 36–52%"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            # Преобразуем в верхний регистр для надежности
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            g_count = seq_upper.count('G')
            c_count = seq_upper.count('C')
            gc_count = g_count + c_count

            if length > 0:
                gc_percent = (gc_count / length) * 100
                g_percent = (g_count / length) * 100
                c_percent = (c_count / length) * 100
            else:
                gc_percent = g_percent = c_percent = 0.0

            gc_point = 0 if gc_percent < 36 else 0 if gc_percent > 52 else 1

            results.append({
                'sequence_id': f"seq_{idx + 1:03d}",
                'sense': sequence,
                'length': length,
                'gc_percent': round(gc_percent, 2),
                'gc_point': gc_point,
            })
        db = pd.DataFrame(results)
        return db

    def frequency_gc_at(self):
        """Frequency of GC and AT repeats less than 3 and 4, respectively"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            gc_frequency = 0
            au_frequency = 0

            # Правильный подсчет со смещением 1
            for i in range(length - 1):  # -1 потому что берем по 2 символа
                if seq_upper[i:i + 2] == 'GC':
                    gc_frequency += 1
                elif seq_upper[i:i + 2] == 'AU':
                    au_frequency += 1

            frequency_point = 0 if gc_frequency > 3 else 0 if au_frequency > 4 else 1
            results.append({
                'sense': sequence,
                'gc_frequency': gc_frequency,
                'au_frequency': au_frequency,
                'frequency_point': frequency_point
            })
        db = pd.DataFrame(results)
        return db

    def add_dt_overhangs(self, dt_count=2):
        """TT overhangs in 3′-end"""
        results = []
        dt_overhang = 'T' * dt_count
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)
            sequence_with_dt = dt_overhang + seq_upper + dt_overhang
            original_length = len(seq_upper)
            extended_length = len(sequence_with_dt)
            dt_point = 1
            results.append({
            'sense_with_dt': sequence_with_dt,
            'original_length': original_length,
            'extended_length': extended_length,
            'dt_point': dt_point
            })
        db = pd.DataFrame(results)
        return db

    def u_at_10(self):
        """Presence of U nucleotide at the tenth position of the sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)
            tenth_nucleotide = seq_upper[9]  # 10-я позиция = индекс 9
            has_u_at_10 = (tenth_nucleotide == 'U')
            u_tenth_point = 1 if has_u_at_10 else 0
            results.append({
                'sense': sequence,
                'u_10_point': u_tenth_point,
            })
        db = pd.DataFrame(results)
        return db

    def no_g_at_13(self):
        """Absence of G nucleotide at the thirteenth position of the sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)
            thirteenth_nucleotide = seq_upper[12]
            hasnt_g_at_13 = (thirteenth_nucleotide != 'G')
            no_g_thirteenth_point = 1 if hasnt_g_at_13 else 0
            results.append({
                'sense': sequence,
                'no_g_13_point': no_g_thirteenth_point,
            })
        db = pd.DataFrame(results)
        return db

    def no_g_c_at_19(self):
        """Absence of G/C at the nineteenth position of the sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            if length >= 19:
                nineteenth_nucleotide = seq_upper[18]
                no_g_c_at_19 = (nineteenth_nucleotide != 'G') and (nineteenth_nucleotide != 'C')
                no_g_c_at_19 = 1 if no_g_c_at_19 else 0
            else:
                no_g_c_at_19 = 0
            results.append({
                'sense': sequence,
                'no_g_c_19_point': no_g_c_at_19,
            })
        db = pd.DataFrame(results)
        return db

    def a_at_3_19(self):
        """Presence of A nucleotide at the third and the nineteenth position of the sense strand"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            if length >= 19:
                nineteenth_nucleotide = seq_upper[18]
                third_nucleotide = seq_upper[2]
                a_at_3_19 = (third_nucleotide == 'A') and (nineteenth_nucleotide == 'A')
                a_at_3_19 = 1 if a_at_3_19 else 0
            else:
                a_at_3_19 = 0
            results.append({
                'sense': sequence,
                'a_at_3_19_point': a_at_3_19,
            })
        db = pd.DataFrame(results)
        return db

    def strong_pairing_5_end(self):
        """Strong base pairing at 5′-end of the sense strand (presence of G/C)"""
        results = []
        for idx, sequence in enumerate(self.df_sense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            if length >= 1:
                first_nucleotide = seq_upper[0]
                g_c_5_end = (first_nucleotide == 'G') or (first_nucleotide == 'C')
                g_c_5_end = 1 if g_c_5_end else 0
            else:
                g_c_5_end = 0
            results.append({
                'sense': sequence,
                'g_c_5_end_point': g_c_5_end,
            })
        db = pd.DataFrame(results)
        return db

    """Обработка antisense-последовательности"""

    def a_at_6(self):
        """Presence of A nucleotide at the sixth position of the antisense strand"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)
            sixth_nucleotide = seq_upper[5]  # 10-я позиция = индекс 9
            a_at_6 = (sixth_nucleotide == 'A')
            a_at_6_point = 1 if a_at_6 else 0
            results.append({
            'antisense': sequence,
            'a_at_6_point': a_at_6_point,
            })
        db = pd.DataFrame(results)
        return db

    def weak_pairing_5_end(self):
        """Weak base pairing at 5′-end of the antisense strand (presence of A/U)"""
        results = []
        for idx, sequence in enumerate(self.df_antisense['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)

            if length >= 1:
                first_nucleotide = seq_upper[0]
                a_u_5_end = (first_nucleotide == 'A') or (first_nucleotide == 'U')
                a_u_5_end = 1 if a_u_5_end else 0
            else:
                a_u_5_end = 0
            results.append({
                'antisense': sequence,
                'a_u_5_end_point': a_u_5_end,
            })
        db = pd.DataFrame(results)
        return db


    """Вывод информации"""

    def final_count(self, columns_to_show=None):
        df1 = self.analyze_gc().add_prefix('m1_')
        df2 = self.frequency_gc_at().add_prefix('m2_')
        df3 = self.add_dt_overhangs().add_prefix('m3_')
        df4 = self.u_at_10().add_prefix('m4_')
        df5 = self.no_g_at_13().add_prefix('m5_')
        df6 = self.no_g_c_at_19().add_prefix('m6_')
        df7 = self.a_at_3_19().add_prefix('m7_')
        df8 = self.strong_pairing_5_end().add_prefix('m8_')
        df9 = self.a_at_6().add_prefix('m9_')
        df10 = self.weak_pairing_5_end().add_prefix('m10_')


        general_df = pd.concat([df1, df2, df3, df4, df5, df6, df7, df8, df9,df10], axis=1)
        point_columns = [col for col in general_df.columns if 'point' in col]
        general_df['total_points'] = general_df[point_columns].sum(axis=1)
        if columns_to_show is None:
            columns_to_show = ['m1_sense'] + ['m9_antisense'] + point_columns + ['total_points']
        else:
            if 'total_points' not in columns_to_show:
                columns_to_show.append('total_points')
            if 'm1_sense' not in columns_to_show:
                columns_to_show = ['m1_sense'] + columns_to_show

        available_cols = [col for col in columns_to_show if col in general_df.columns]
        return general_df[available_cols]

    def get_combined_data(self, columns_to_show=None):
        """Объединяет все три метода и позволяет выбрать столбцы"""
        # Получаем DataFrame от всех методов
        df1 = self.analyze_gc().add_prefix('m1_')
        df2 = self.frequency_gc_at().add_prefix('m2_')
        df3 = self.add_dt_overhangs().add_prefix('m3_')

        combined_df = pd.concat([df1, df2, df3], axis=1)

        if columns_to_show:
            available_columns = [col for col in columns_to_show if col in combined_df.columns]
            combined_df = combined_df[available_columns]
        return combined_df



d = analyz_rna()
# g = d.analyze_gc()
# c = helper()
# print(g.to_csv("my_data.csv"))
# h = d.a_at_6()
# print(h.to_csv("a_at_6.csv"))
# h = d.weak_pairing_5_end()
# print(h.to_csv("weak_pairing_5_end.csv"))
# t = d.get_combined_data([])
t = d.final_count()
print(t.to_csv('general_data.csv'))