import pandas as pd
from prepare_rna import df_sense


def convert_chromosome_format(chrom):
    if pd.isna(chrom): return chrom
    chrom_str = str(chrom)
    if chrom_str.startswith('NC_0000'):
        if 'chr' in chrom_str:
            parts = chrom_str.split(':')
            if len(parts) > 1 and 'chr' in parts[1]:
                return parts[1].replace('chr', '')
        else:
            return chrom_str.split('.')[0][-2:].lstrip('0')
    return chrom_str.replace('chr', '')


def read_bed_file(file_path):
    bed_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    df = pd.read_csv(file_path, sep='\t', header=None, comment='#')
    df.columns = bed_columns[:len(df.columns)]
    df['chrom'] = df['chrom'].apply(convert_chromosome_format)
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end'] = pd.to_numeric(df['end'], errors='coerce')
    return df


def is_critical_position(pos, sirna_length):
    # Критические позиции адаптированы под длину siRNA
    critical_positions = {
        1: "5'-конец",
        19: "3'-конец",
        10: "сайт разрезания",
        11: "сайт разрезания"
    }
    # Seed регион (позиции 2-8) для любой длины
    if 2 <= pos <= 8:
        return True, "seed-регион"
    return pos in critical_positions, critical_positions.get(pos, "не критичная")


def analyze_sirna_snp(sequence, start_pos, snps_in_exon, exon_start, sirna_length):
    result = {'has_snp': False, 'total_snps': 0, 'critical_snps': 0, 'snp_names': []}

    for _, snp in snps_in_exon.iterrows():
        snp_pos_in_exon = snp['start'] - exon_start
        if start_pos <= snp_pos_in_exon < start_pos + sirna_length:
            result['has_snp'] = True
            result['total_snps'] += 1
            snp_pos_in_sirna = snp_pos_in_exon - start_pos
            is_critical, _ = is_critical_position(snp_pos_in_sirna, sirna_length)
            if is_critical:
                result['critical_snps'] += 1
                result['snp_names'].append(snp['name'])

    return result


def create_sirna_dataframe(sequence_df, snp_df, exon_coords):
    snps_in_exon = snp_df[
        (snp_df['chrom'] == str(exon_coords['chrom'])) &
        (snp_df['start'] >= exon_coords['start']) &
        (snp_df['end'] <= exon_coords['end'])
        ]

    data = []
    for _, seq_row in sequence_df.iterrows():
        sequence = str(seq_row['sequence']).upper().replace('U', 'T')
        sirna_length = seq_row['size_nt']  # Используем реальную длину из данных

        # Генерируем siRNA кандидаты скользящим окном
        for i in range(len(sequence) - sirna_length + 1):
            snp_analysis = analyze_sirna_snp(sequence, i, snps_in_exon, exon_coords['start'], sirna_length)

            not_in_snp_score = 1 if not snp_analysis['has_snp'] else 0
            no_critical_score = 1 if snp_analysis['critical_snps'] == 0 else 0
            snp_avoidance_score = 1 if not snp_analysis['has_snp'] else 0

            data.append({
                'fragment_id': seq_row['fragment_id'],
                'sirna_sequence': sequence[i:i + sirna_length],
                'sirna_length': sirna_length,
                'has_snp': snp_analysis['has_snp'],
                'total_snps': snp_analysis['total_snps'],
                'critical_snps': snp_analysis['critical_snps'],
                'not_in_snp_sites_score': not_in_snp_score,
                'no_critical_snps_score': no_critical_score,
                'snp_avoidance_score': snp_avoidance_score,
                'total_score': snp_avoidance_score,
                'snp_names': ', '.join(snp_analysis['snp_names'])
            })

    return pd.DataFrame(data)


# Основной код
snp_df = read_bed_file("Live RefSNPs dbSNP b157 v2.BED")
exon_coords = {'chrom': '6', 'start': 16326394, 'end': 16328470}

results_df = create_sirna_dataframe(df_sense, snp_df, exon_coords)
results_df.to_csv('sirna_snp_results.csv', index=False)
print("Анализ SNP\n"
'has_snp'        " - "   "Есть ли хотя бы один SNP в этой siRNA (True/False)\n"
'total_snps'     " - "    "Общее количество SNP в siRNA\n"
'critical_snps'  " - "    "Количество SNP в критических позициях\n"
'snp_names'      " - "    "Имена критических SNP через запятую\n"

"Система баллов\n"
'not_in_snp_sites_score'  " - "   "1 если нет SNP, 0 если есть SNP\n"
'no_critical_snps_score'  " - "   "1 если нет критических SNP, 0 если есть\n"
'snp_avoidance_score'     " - "   "1 (нет SNP), 0 (есть критические SNP, есть SNP но не критические)\n"
'total_score'             " - "   "Общий балл (равен snp_avoidance_score)\n"
)