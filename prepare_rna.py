from datetime import datetime
import pandas as pd
from Bio import SeqIO
sequence1 = SeqIO.read("sequence.fasta", "fasta")
ATXN1 = str(sequence1.seq)


class editor_rna:
    """Метод класса для работы с размерами RNA"""
    def __init__(self):
        self.sequence = ATXN1

    def info_rna(self):
        """Метод выводит последовательность РНК"""
        print(self.sequence)
        print(f"Длина последовательности: {len(self.sequence)}")

    def sense(self, start_size=15, end_size=30, filename=None):
        sense = self.sequence.replace('T', 'U')
        sense = sense[781:2858]
        start_size = 15
        end_size = 30 + 1
        fragments_list = []
        print(f"Длина фрагмента для нарезки: {len(sense)} нуклеотидов")
        for fragment_size in range(start_size, end_size):
            for start_pos in range(0, len(sense), fragment_size):
                end_pos = start_pos + fragment_size
                if end_pos > len(sense):
                    break
                fragment = sense[start_pos:end_pos]

                fragments_list.append({
                    'fragment_id': f"{fragment_size}_{start_pos + 1}",
                    'size_nt': fragment_size,
                    'sequence': fragment,
                    'sequence_length': len(fragment)
                })
        df_sense = pd.DataFrame(fragments_list)
        if filename:
            df_sense.to_csv(filename, index=False)
            print(f"Результаты сохранены в файл: {filename}")
        print(f"Создано {len(fragments_list)} фрагментов")
        return df_sense

    def antisense(self, start_size=15, end_size=30, filename=None):
        sense = self.sequence.replace('T', 'U')
        comp_dict = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
        antisense = ''.join(comp_dict.get(nuc, nuc) for nuc in sense.upper())
        antisense = antisense[781:2858]
        start_size = 15
        end_size = 30 + 1
        fragments_list = []
        print(f"Длина фрагмента для нарезки: {len(antisense)} нуклеотидов")
        for fragment_size in range(start_size, end_size):
            for start_pos in range(0, len(antisense), fragment_size):
                end_pos = start_pos + fragment_size
                if end_pos > len(antisense):
                    break
                fragment = antisense[start_pos:end_pos]

                fragments_list.append({
                    'fragment_id': f"{fragment_size}_{start_pos + 1}",
                    'size_nt': fragment_size,
                    'sequence': fragment,
                    'sequence_length': len(fragment)
                })
        df_antisense = pd.DataFrame(fragments_list)
        if filename:
            df_antisense.to_csv(filename, index=False)
            print(f"Результаты сохранены в файл: {filename}")
        print(f"Создано {len(fragments_list)} фрагментов")
        return df_antisense

    def slice_rna(self, start_size=15, end_size=30, filename=None):
        # start = int(input("Введите начальный индекс: "))
        # end = int(input("Введите конечный индекс: "))
        # frag = self.sequence[start:end]
        frag = self.sequence[781:2858]
        # start_size = int(input("Введите начальный размер: "))
        # end_size = int(input("Введите конечный размер: ")) + 1
        start_size = 15
        end_size = 30 + 1
        fragments_list = []
        print(f"Длина фрагмента для нарезки: {len(frag)} нуклеотидов")
        for fragment_size in range(start_size, end_size):
            for start_pos in range(0, len(frag), fragment_size):
                end_pos = start_pos + fragment_size
                if end_pos > len(frag):
                    break
                fragment = frag[start_pos:end_pos]

                fragments_list.append({
                    'fragment_id': f"{fragment_size}_{start_pos + 1}",
                    'size_nt': fragment_size,
                    'sequence': fragment,
                    'sequence_length': len(fragment)
                })
        df = pd.DataFrame(fragments_list)
        if filename:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            filename = f"analyz_rna_{timestamp}.csv"
            df.to_csv(filename, index=False)
            print(f"Результаты сохранены в файл: {filename}")
        print(f"Создано {len(fragments_list)} фрагментов")
        return df

a = editor_rna()
df_sense = a.sense()
df_antisense = a.antisense()


