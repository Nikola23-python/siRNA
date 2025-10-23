def u_at_tenth(sequence: str):
    """
    Check if U nucleotide is at the tenth position (1-based indexing).

    Args:
        sequence: RNA sequence string

    Returns:
        bool: True if U at position 10, False otherwise
    """
    if sequence[9].upper() == 'U':
        return  1
    return 0


# Usage
rna_sequence = "AUGCCGACGGUACGUAG"
result = has_u_at_tenth(rna_sequence)
print(f"U at tenth position: {result}")

# from analyze_rna import analyze_rna
#
# #[781:2858]
# #info_about_methods - поможет понять, что есть в классе analyze_rna
#
# from Bio import SeqIO
#
# sequence = SeqIO.read("sequence.fasta", "fasta")
# seq = sequence.seq
# siRNA = analyze_rna(sequence)
# seq = siRNA.range_sequence()
#
#
# import pandas as pd
# from datetime import datetime
#
#
# def slice_rna_to_csv(sequence, start_size=15, end_size=30, filename=None):
#
#     if filename is None:
#         timestamp = datetime.now().strftime("%Y%m%d_%H%M")
#         filename = f"rna_fragments_{timestamp}.csv"
#
#     fragments_list = []
#
#     print("🔪 Нарезка РНК последовательности...")
#
#     for fragment_size in range(start_size, end_size):
#         for start_pos in range(0, len(sequence), fragment_size):
#             end_pos = start_pos + fragment_size
#             if end_pos > len(sequence):
#                 break
#             fragment = sequence[start_pos:start_pos + fragment_size]
#             nucleotide_counts = {
#                 'A': fragment.count('A'),
#                 'U': fragment.count('U'),
#                 'G': fragment.count('G'),
#                 'C': fragment.count('C')
#             }
#             fragments_list.append({
#                 'fragment_id': f"{fragment_size}_{start_pos + 1}",
#                 'size_nt': fragment_size,
#                 'sequence': fragment,
#                 **nucleotide_counts,  # Распаковываем счетчики нуклеотидов
#                 'GC_content': round((nucleotide_counts['G'] + nucleotide_counts['C']) / fragment_size * 100, 2),
#                 'sequence_length': len(sequence)
#             })
#
#     df = pd.DataFrame(fragments_list)
#
#     df.to_csv(filename, index=False, encoding='utf-8')
#
#     print(f"✅ Данные сохранены в: {filename}")
#     print(f"📊 Статистика:")
#     print(f"   - Всего фрагментов: {len(df)}")
#     print(f"   - Диапазон размеров: {start_size}-{end_size} нт")
#     print(f"   - Размеры в файле: {df['size_nt'].unique()}")
#
#     return df
#
# from Bio import SeqIO
#
# sequence = SeqIO.read("sequence.fasta", "fasta")
# seq = sequence.seq
# siRNA = analyze_rna(sequence)
# seq1 = siRNA.range_sequence()
# print(type(seq1))
# slice_rna_to_csv(seq1)