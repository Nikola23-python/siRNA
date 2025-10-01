from datetime import datetime

import pandas as pd
from editor_rna import df
from helper import helper


class analyz_rna:
    def __init__(self):
        self.df = df

    def analyze_gc(self, filename=None):
        """GC content of 36–52%"""
        results = []

        for idx, sequence in enumerate(self.df['sequence']):
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
                'sequence': sequence,
                'length': length,
                # 'g_count': g_count,
                # 'c_count': c_count,
                # 'gc_count': gc_count,
                # 'g_percent': round(g_percent, 2),
                # 'c_percent': round(c_percent, 2),
                'gc_percent': round(gc_percent, 2),
                'gc_point': gc_point,
            })
        db = pd.DataFrame(results)
        return db

    def frequency_gc_at(self, filename=None):
        """Frequency of GC and AT repeats less than 3 and 4, respectively"""
        b = analyz_rna()
        a = b.analyze_gc()
        for idx, sequence in enumerate(self.df['sequence']):
            seq_upper = str(sequence).upper()
            length = len(seq_upper)
            gc_frequency = seq_upper.count('GC')
            a['gc_frequency'] = gc_frequency
            at_frequency = seq_upper.count('AT')
            a['at_frequency'] = at_frequency
            a['frequency_point'] = 0 if gc_frequency > 3 else 0 if at_frequency > 4 else 1
            # if filename is None:
            #     timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            #     filename = f"analyz_rna_{timestamp}.csv"
            # a.to_csv(filename, index=False, encoding='utf-8')
            return a

    def add_dt_overhangs(self, dt_count=2, filename=None):
        """TT overhangs in 3′-end"""
        b = analyz_rna()
        a = b.frequency_gc_at()
        dt_overhang = 'T' * dt_count
        a['sequence_with_dt'] = a['sequence'].apply(
                lambda x: dt_overhang + str(x) + dt_overhang)
        a['original_length'] = a['sequence'].str.len()
        a['extended_length'] = a['sequence_with_dt'].str.len()
        a['dt_point'] = 1
        print(f"✅ Добавлено {dt_count} dT с каждой стороны к {len(a)} последовательностям")
        if filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            filename = f"analyz_rna_{timestamp}.csv"
        a.to_csv(filename, index=False, encoding='utf-8')
        return a


d = analyz_rna()
c = d.add_dt_overhangs()