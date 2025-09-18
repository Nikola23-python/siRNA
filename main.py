from Bio import SeqIO
from Bio.Seq import Seq

sequence = SeqIO.read("sequence.fasta", "fasta")
seq = sequence.seq

seq1 = seq[781:2858]
print(seq1)

def split_rna_19(sequence):
    fragment_length = 19
    fragments = [sequence[i:i + fragment_length] for i in range(0, len(sequence), fragment_length)]
    return fragments

seq2 = list(split_rna_19(seq1))
print(seq2)

# # print(f"Всего фрагментов: {len(result)}")
# # for i, fragment in enumerate(result, 1):
# #     print(f"Фрагмент {i}: {fragment} (длина: {len(fragment)})")
#
#
def calculate_gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    total_length = len(sequence)
    if total_length > 0:
        gc_percent = (gc_count / total_length) * 100
    else:
        gc_percent = 0.0

    return gc_percent

sequence = 'CCTGGATCCCCTACAGAAA'
print(calculate_gc_content(sequence))

# gc_percentages = [calculate_gc_content(seq) for seq in seq2]
#
# # Выводим результаты
# for i, (seq, gc_percent) in enumerate(zip(seq2, gc_percentages)):
#     print(f"Последовательность {i+1}: '{seq[:15]}...' | GC: {gc_percent:.1f}%")

import pandas as pd


def analyze_sequences_gc(seq2):
    """Анализирует GC-content для списка последовательностей и возвращает DataFrame"""
    results = []

    for i, seq in enumerate(seq2):
        gc_percent = calculate_gc_content(seq)
        length = len(seq)
        gc_count = seq.upper().count('G') + seq.upper().count('C')
        at_count = length - gc_count

        results.append({
            'sequence_id': f'seq_{i + 1:03d}',
            'sequence_preview': seq[:20] + '...' if len(seq) > 20 else seq,
            'length': length,
            'gc_count': gc_count,
            'at_count': at_count,
            'gc_percent': gc_percent,
            'gc_category': 'Low' if gc_percent < 36 else 'High' if gc_percent > 52 else 'Optimal'
        })

    return pd.DataFrame(results)



# Создаем DataFrame с результатами
df = analyze_sequences_gc(seq2)

# Выводим красивые результаты
print("Результаты анализа GC-content:")
print(df[['sequence_id', 'length', 'gc_percent', 'gc_category']].to_string(index=False))

# Статистика
print(f"\nСтатистика:")
print(f"Средний GC-content: {df['gc_percent'].mean():.1f}%")
print(f"Минимальный GC-content: {df['gc_percent'].min():.1f}%")
print(f"Максимальный GC-content: {df['gc_percent'].max():.1f}%")

with pd.ExcelWriter('gc_analys.xlsx') as writer:
    # Основные данные
    df.to_excel(writer, sheet_name='Sequences', index=False)

    # Статистика
    df.to_excel(writer, sheet_name='Statistics', index=False)

print("Файл 'gc_analys.xlsx' успешно создан!")