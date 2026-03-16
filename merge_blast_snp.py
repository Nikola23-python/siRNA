import pandas as pd
import os
import glob
from datetime import datetime

def load_blast_results():
    pattern = 'BLAST/BLAST_analysis.csv'
    blast_df = pd.read_csv(pattern)
    return blast_df

def load_snp_results():
    pattern = 'SNP/SNP_analysis_.csv'
    snp_df = pd.read_csv(pattern)
    return snp_df

def prepare_data(blast_df, snp_df):
    blast_df['fragment_id'] = blast_df['fragment_id'].astype(str)
    snp_df['fragment_id'] = snp_df['fragment_id'].astype(str)
    blast_ids = set(blast_df['fragment_id'])
    snp_ids = set(snp_df['fragment_id'])
    common_ids = blast_ids.intersection(snp_ids)
    if len(common_ids) == 0:
        print("\n⚠️  ВНИМАНИЕ: Нет общих fragment_id!")
        print("   Примеры из BLAST:", list(blast_ids)[:5])
        print("   Примеры из SNP:", list(snp_ids)[:5])
        return None, None, None
    return blast_ids, snp_ids, common_ids

def merge_results(blast_df, snp_df):
    blast_ids, snp_ids, common_ids = prepare_data(blast_df, snp_df)
    if common_ids is None:
        return None
    snp_columns = ['fragment_id', 'has_snp', 'snp_score', 'num_snps']
    optional_cols = ['top_snp', 'top_snp_maf', 'genomic_start', 'genomic_end', 'size_nt']
    for col in optional_cols:
        if col in snp_df.columns:
            snp_columns.append(col)
    merged = pd.merge(
        blast_df,
        snp_df[snp_columns],
        on='fragment_id',
        how='inner',
        suffixes=('', '_snp')
    )
    if 'length' not in merged.columns:
        merged['length'] = merged['size_nt']
    if 'antisense_sequence' in merged.columns:
        merged['antisense_sequence'] = merged['antisense_sequence'].str.replace('T', 'U')
    merged['total_score'] = merged['blast_score'] + merged['snp_score']

    def get_category(row):
        blast = row['blast_score']
        snp = row['snp_score']
        if blast == 2 and snp == 1:
            return 'Excellent'
        elif blast == 2 and snp == 0:
            return 'Good (SNP inside)'
        elif blast == 1 and snp == 1:
            return 'Good (Off-target risk)'
        else:
            return 'Poor'
    merged['category'] = merged.apply(get_category, axis=1)
    return merged

def save_merged_results(merged_df):
    all_file = 'sirna_merged_all.csv'
    merged = merged_df[merged_df['category'] == 'Excellent']
    merged.to_csv(all_file, index=False)
    print(f"   ✓ Все siRNA: {all_file}")
    print(f"     Записей: {len(merged_df)}")
    return all_file

def print_summary(merged_df):
    total = len(merged_df)
    print(f"\n💡 РЕКОМЕНДАЦИИ:")
    if 'category' in merged_df.columns:
        excellent_count = len(merged_df[merged_df['category'] == 'Excellent'])
        good_count = len(merged_df[merged_df['category'] == 'Good'])
        print(f"   • Отличных кандидатов (Excellent): {excellent_count}")
        print(f"   • Хороших кандидатов (Good): {good_count}")
        print(
            f"   • Всего пригодных: {excellent_count + good_count} ({(excellent_count + good_count) / total * 100:.1f}%)")

def main():
    blast_df = load_blast_results()
    snp_df = load_snp_results()
    merged_df = merge_results(blast_df, snp_df)

    if merged_df is None or len(merged_df) == 0:
        print("\n❌ Ошибка: нет данных для объединения")
        return

    save_merged_results(merged_df)
    print_summary(merged_df)

if __name__ == "__main__":
    main()