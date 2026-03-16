import pandas as pd
import subprocess
import os
import multiprocessing as mp
from collections import defaultdict
from tqdm import tqdm
import time

BLAST_DB = "/home/nikolay/PycharmProjects/mRNA/BLAST/blast_dbs2/human_refseq"
BATCH_SIZE = 500
THREADS_PER_BLAST = 2
NUM_WORKERS = mp.cpu_count() // 2

def run_single_batch(batch_data):
    batch_idx, sequences = batch_data
    tmp_file = f"temp_batch_{batch_idx}.fasta"
    with open(tmp_file, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq_{i}\n{seq}\n")
    cmd = [
        "blastn",
        "-query", tmp_file,
        "-db", BLAST_DB,
        "-task", "blastn-short",
        "-word_size", "7",
        "-evalue", "3000",
        "-perc_identity", "78",
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle qcovs",
        "-num_threads", str(THREADS_PER_BLAST),
        "-dust", "no"
    ]
    results = defaultdict(list)
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if res.stdout.strip():
            for line in res.stdout.strip().split('\n'):
                parts = line.split('\t')
                results[parts[0]].append(line)
        status = "success"
    except subprocess.TimeoutExpired:
        status = "timeout"
    except Exception as e:
        status = f"error: {str(e)}"
    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
    return results, sequences, status

def analyze_specificity(hits, sequence, target="ATXN1"):
    if hits is None:
        return False, "BLAST Error/Timeout"
    if not hits:
        return True, "Perfect match to ATXN1"
    seq_len = len(sequence)
    for line in hits:
        parts = line.split('\t')
        title = parts[12].upper()
        pident = float(parts[2])
        aln_len = int(parts[3])
        if target.upper() not in title:
            if pident >= 78 and aln_len >= 15:
                off_target_name = parts[12][:50]
                return False, f"Off-target: {off_target_name}"
    return True, "Perfect match to ATXN1"


def main_parallel_blast(df_sirna):
    unique_seqs = pd.concat([df_sirna['sense_sequence'], df_sirna['antisense_sequence']]).unique().tolist()
    batches = [(i, unique_seqs[i:i + BATCH_SIZE]) for i in range(0, len(unique_seqs), BATCH_SIZE)]
    print(f"🚀 Запуск параллельного BLAST: {len(batches)} батчей на {NUM_WORKERS} воркерах")
    all_seq_results = {}
    with mp.Pool(NUM_WORKERS) as pool:
        for batch_hits, batch_seqs, status in tqdm(pool.imap_unordered(run_single_batch, batches), total=len(batches)):
            for i, seq in enumerate(batch_seqs):
                if status == "success":
                    hits = batch_hits.get(f"seq_{i}", [])
                    # analyze_specificity теперь возвращает (True/False, "Reason")
                    all_seq_results[seq] = analyze_specificity(hits, seq)
                else:
                    all_seq_results[seq] = (False, f"Batch Error: {status}")
    print("📊 Финализация данных...")
    df_sirna['sense_ok'] = df_sirna['sense_sequence'].map(lambda x: all_seq_results[x][0])
    df_sirna['sense_reason'] = df_sirna['sense_sequence'].map(lambda x: all_seq_results[x][1])
    df_sirna['anti_ok'] = df_sirna['antisense_sequence'].map(lambda x: all_seq_results[x][0])
    df_sirna['anti_reason'] = df_sirna['antisense_sequence'].map(lambda x: all_seq_results[x][1])
    df_sirna['blast_score'] = df_sirna['sense_ok'].astype(int) + df_sirna['anti_ok'].astype(int)
    os.makedirs("BLAST", exist_ok=True)
    output_file = "BLAST/BLAST_analysis.csv"
    df_sirna.to_csv(output_file, index=False)
    print(f"\nРезультаты сохранены в файл: {output_file}")
    print(f"   Специфичных дуплексов (Score 2): {len(df_sirna[df_sirna['blast_score'] == 2])}")
    print(f"   Частично специфичных (Score 1): {len(df_sirna[df_sirna['blast_score'] == 1])}")
    print(f"   Неспецифичных / Ошибки (Score 0): {len(df_sirna[df_sirna['blast_score'] == 0])}")
    return df_sirna

if __name__ == "__main__":
    from prepare_rna import siRNAEditor, ATXN1
    editor = siRNAEditor(ATXN1)
    df_pairs = editor.generate_pairs()
    final_df = main_parallel_blast(df_pairs)

