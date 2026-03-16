from Bio import SeqIO
sequence1 = SeqIO.read("/home/nikolay/PycharmProjects/mRNA/ATXN1/sequence.fasta", "fasta")
ATXN1 = str(sequence1.seq)
import pandas as pd
from Bio.Seq import Seq


class siRNAEditor:
    def __init__(self, sequence):
        self.sequence = sequence.upper().replace('T', 'U')

    def generate_pairs(self, start_idx=781, end_idx=2858, min_len=19, max_len=25):
        target_region = self.sequence[start_idx:end_idx]
        results = []

        for length in range(min_len, max_len + 1):
            for i in range(len(target_region) - length + 1):
                sense_seq = target_region[i: i + length]
                anti_seq = str(Seq(sense_seq).reverse_complement())
                results.append({
                    'fragment_id': f"{length}_{i + 1}",
                    'size_nt': length,
                    'sense_sequence': sense_seq,
                    'antisense_sequence': anti_seq,
                    'genomic_start': 16326394 + start_idx + i,
                    'genomic_end': 16326394 + start_idx + i + length
                })
        df = pd.DataFrame(results)
        print(f"Создано пар siRNA: {len(df)}")
        return df

editor = siRNAEditor(ATXN1)
df_final = editor.generate_pairs()