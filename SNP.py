import pandas as pd
from prepare_rna import siRNAEditor, ATXN1

def analyze_snps_for_sirna(df_sirna):
    EXON_START = 16326394
    MIN_MAF = 0.01
    gnomad = pd.read_csv("SNP/gnomAD_v4.1.0_ENSG00000124788_2025_12_18_01_08_51.csv")
    gnomad['Position'] = pd.to_numeric(gnomad['Position'], errors='coerce')
    gnomad['Allele Frequency'] = pd.to_numeric(gnomad['Allele Frequency'], errors='coerce')
    frequent_snps = gnomad[gnomad['Allele Frequency'] >= MIN_MAF].copy()
    snp_map = frequent_snps.set_index('Position')['Allele Frequency'].to_dict()
    rsid_map = frequent_snps.set_index('Position')['rsIDs'].to_dict()

    def check_row(row):
        start = row['genomic_start']
        end = row['genomic_end']
        found_snps = []
        for pos in range(start, end + 1):
            if pos in snp_map:
                found_snps.append(f"{rsid_map[pos]}({snp_map[pos] * 100:.1f}%)")
        return pd.Series({
            'has_snp': len(found_snps) > 0,
            'snp_score': 1 if len(found_snps) == 0 else 0,
            'num_snps': len(found_snps),
            'top_snp': "; ".join(found_snps) if found_snps else "None"
        })
    snp_results = df_sirna.apply(check_row, axis=1)
    final_df = pd.concat([df_sirna, snp_results], axis=1)
    return final_df

if __name__ == "__main__":
    editor = siRNAEditor(ATXN1)
    df_pairs = editor.generate_pairs()
    final_results = analyze_snps_for_sirna(df_pairs)
    output_path = "SNP/SNP_analysis_.csv"
    final_results.to_csv(output_path, index=False)

    print(f"\nСохранено в: {output_path}")
    print(f"Чистых кандидатов (Score 1): {len(final_results[final_results['snp_score'] == 1])}")