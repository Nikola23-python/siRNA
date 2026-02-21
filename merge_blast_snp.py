import pandas as pd
import os
import glob
from pathlib import Path

def find_latest_snp_file():
    snp_files = glob.glob('SNP/all_fragments_*_with_snps_*.csv')
    if not snp_files:
        raise FileNotFoundError("Не найдены файлы SNP анализа в папке SNP/")
    # Сортируем по времени модификации (самый свежий последний)
    latest_file = max(snp_files, key=os.path.getctime)
    return latest_file

def load_blast_results(blast_file='sirna_blast_good_results.csv'):
    if not os.path.exists(blast_file):
        blast_files = glob.glob('BLAST/sirna_blast_*results.csv')
        if not blast_files:
            raise FileNotFoundError("Не найдены файлы BLAST результатов")
        blast_file = max(blast_files, key=os.path.getctime)
    print(f"📊 Загрузка BLAST результатов из: {blast_file}")
    blast_df = pd.read_csv(blast_file)
    print(f"   Загружено {len(blast_df)} записей")
    return blast_df

def load_snp_results(snp_file=None):
    if snp_file is None:
        snp_file = find_latest_snp_file()

    print(f"🧬 Загрузка SNP результатов из: {snp_file}")
    snp_df = pd.read_csv(snp_file)
    print(f"   Загружено {len(snp_df)} записей")
    return snp_df


def prepare_snp_data(snp_df):
    snp_prepared = snp_df.copy()
    snp_prepared['base_id'] = snp_prepared['fragment_id']
    snp_prepared['is_sirna_length'] = snp_prepared['length'].between(19, 23)
    snp_prepared['snp_score'] = 0
    snp_prepared.loc[~snp_prepared['has_snp'], 'snp_score'] = 2  # Нет SNP
    snp_prepared.loc[
        snp_prepared['has_snp'] & ~snp_prepared['critical_snp'], 'snp_score'] = 1  # Есть SNP, но не критические
    return snp_prepared


def prepare_blast_data(blast_df):
    blast_prepared = blast_df.copy()
    blast_prepared['base_id'] = blast_prepared['fragment_id']
    return blast_prepared

def merge_results(blast_df, snp_df):
    print("\n" + "=" * 80)
    print(" ОБЪЕДИНЕНИЕ РЕЗУЛЬТАТОВ BLAST И SNP АНАЛИЗА")
    print("=" * 80)

    blast_prep = prepare_blast_data(blast_df)
    snp_prep = prepare_snp_data(snp_df)

    print(f"\n📊 Исходные данные:")
    print(f"   BLAST: {len(blast_prep)} записей")
    print(f"   SNP: {len(snp_prep)} записей")

    # Находим пересекающиеся fragment_id
    blast_ids = set(blast_prep['base_id'])
    snp_ids = set(snp_prep['base_id'])
    common_ids = blast_ids.intersection(snp_ids)

    print(f"\n🔄 Пересечение:")
    print(f"   В BLAST: {len(blast_ids)} уникальных ID")
    print(f"   В SNP: {len(snp_ids)} уникальных ID")
    print(f"   Общих: {len(common_ids)} ID")

    if len(common_ids) == 0:
        print("⚠ ВНИМАНИЕ: Нет общих fragment_id!")
        print("   Возможно разные форматы идентификаторов")
        print("\nПримеры из BLAST:", list(blast_ids)[:5])
        print("Примеры из SNP:", list(snp_ids)[:5])
        return None

    # Объединение данных
    merged = pd.merge(
        blast_prep,
        snp_prep[['base_id', 'sequence', 'length', 'has_snp', 'critical_snp',
                  'num_snps', 'top_snp', 'top_snp_maf', 'snp_score',
                  'genomic_start', 'genomic_end', 'start_in_exon', 'end_in_exon']],
        on='base_id',
        how='inner',
        suffixes=('_blast', '_snp')
    )

    # Проверяем совпадение последовательностей
    seq_match = (merged['sense_sequence'] == merged['sequence']) | \
                (merged['antisense_sequence'] == merged['sequence'])

    if not seq_match.all():
        mismatch_count = (~seq_match).sum()
        print(f"\n⚠ ВНИМАНИЕ: {mismatch_count} записей имеют несовпадающие последовательности")

        # Для несовпадающих проверяем, может быть это antisense?
        merged['seq_match_sense'] = merged['sense_sequence'] == merged['sequence']
        merged['seq_match_antisense'] = merged['antisense_sequence'] == merged['sequence']

    # Создаем итоговый score
    merged['total_score'] = merged['blast_score'] + merged['snp_score']

    # Добавляем категорию
    def get_category(row):
        if row['blast_score'] == 2 and row['snp_score'] == 2:
            return 'Отличная'
        elif row['blast_score'] == 2 and row['snp_score'] == 1:
            return 'Хорошая (есть SNP)'
        elif row['blast_score'] == 1 and row['snp_score'] == 2:
            return 'Хорошая (одна цепь неспецифична)'
        elif row['blast_score'] == 2 and row['snp_score'] == 0:
            return 'Средняя (критические SNP)'
        elif row['blast_score'] == 0 and row['snp_score'] == 2:
            return 'Средняя (неспецифична)'
        else:
            return 'Плохая'

    merged['category'] = merged.apply(get_category, axis=1)

    return merged


def analyze_merged_results(merged_df):
    """Анализ объединенных результатов"""
    print("\n" + "=" * 80)
    print(" АНАЛИЗ ОБЪЕДИНЕННЫХ РЕЗУЛЬТАТОВ")
    print("=" * 80)

    total = len(merged_df)

    print(f"\n📊 ОБЩАЯ СТАТИСТИКА:")
    print(f"   Всего siRNA после фильтрации: {total}")

    # Статистика по BLAST
    blast_stats = merged_df['blast_score'].value_counts().sort_index()
    print(f"\n🎯 BLAST SCORE:")
    for score, count in blast_stats.items():
        print(f"   Score {score}: {count} ({count / total * 100:.1f}%)")

    # Статистика по SNP
    snp_stats = merged_df['snp_score'].value_counts().sort_index()
    print(f"\n🧬 SNP SCORE:")
    for score, count in snp_stats.items():
        if score == 2:
            desc = "нет SNP"
        elif score == 1:
            desc = "некритические SNP"
        else:
            desc = "критические SNP"
        print(f"   Score {score} ({desc}): {count} ({count / total * 100:.1f}%)")

    # Распределение по категориям
    print(f"\n🏷️ КАТЕГОРИИ:")
    category_stats = merged_df['category'].value_counts()
    for category, count in category_stats.items():
        print(f"   {category}: {count} ({count / total * 100:.1f}%)")

    # Распределение по длинам
    print(f"\n📏 РАСПРЕДЕЛЕНИЕ ПО ДЛИНАМ:")
    length_stats = merged_df['length'].value_counts().sort_index()
    for length, count in length_stats.items():
        print(f"   {length} нт: {count} ({count / total * 100:.1f}%)")

    # Топ лучших siRNA
    print(f"\n🏆 ТОП-10 ЛУЧШИХ siRNA (max total_score):")
    best_sirnas = merged_df.nlargest(10, 'total_score')
    for idx, row in best_sirnas.iterrows():
        print(f"\n   {row['fragment_id']} ({row['length']}нт):")
        print(f"     Sense: {row['sense_sequence']}")
        print(f"     Anti:  {row['antisense_sequence']}")
        print(f"     BLAST: {row['blast_score']}, SNP: {row['snp_score']}, Total: {row['total_score']}")
        if pd.notna(row.get('top_snp')):
            print(f"     SNP: {row['top_snp']}")

    return category_stats


def save_merged_results(merged_df, prefix='merged_sirna'):
    """Сохранение объединенных результатов"""
    timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")

    # Все результаты
    all_file = f'{prefix}_all_{timestamp}.csv'
    merged_df.to_csv(all_file, index=False)
    print(f"\n💾 Все результаты сохранены в: {all_file}")

    # Только лучшие (total_score >= 3)
    best_file = f'{prefix}_best_{timestamp}.csv'
    best_df = merged_df[merged_df['total_score'] >= 3]
    best_df.to_csv(best_file, index=False)
    print(f"💾 Лучшие siRNA (score >=3) сохранены в: {best_file}")
    print(f"   Количество: {len(best_df)}")

    # Отличные siRNA (total_score == 4)
    excellent_file = f'{prefix}_excellent_{timestamp}.csv'
    excellent_df = merged_df[merged_df['total_score'] == 4]
    excellent_df.to_csv(excellent_file, index=False)
    print(f"💾 Отличные siRNA (score=4) сохранены в: {excellent_file}")
    print(f"   Количество: {len(excellent_df)}")

    # Сводная статистика
    stats = merged_df.groupby(['blast_score', 'snp_score']).size().reset_index(name='count')
    stats_file = f'{prefix}_stats_{timestamp}.csv'
    stats.to_csv(stats_file, index=False)
    print(f"💾 Статистика сохранена в: {stats_file}")

    return all_file, best_file, excellent_file


def main():
    """Основная функция"""
    print("=" * 80)
    print(" ОБЪЕДИНЕНИЕ BLAST И SNP РЕЗУЛЬТАТОВ")
    print("=" * 80)

    try:
        # Загрузка данных
        blast_df = load_blast_results()
        snp_df = load_snp_results()

        # Объединение
        merged_df = merge_results(blast_df, snp_df)

        if merged_df is None or len(merged_df) == 0:
            print("❌ Ошибка: нет данных для объединения")
            return

        # Анализ
        analyze_merged_results(merged_df)

        # Сохранение
        save_merged_results(merged_df)

        print("\n✅ ОБЪЕДИНЕНИЕ ЗАВЕРШЕНО УСПЕШНО!")

    except FileNotFoundError as e:
        print(f"❌ Ошибка: {e}")
        print("\nПроверьте наличие файлов:")
        print("  - BLAST результаты (sirna_blast_*results.csv)")
        print("  - SNP результаты (SNP/all_fragments_*_with_snps_*.csv)")
    except Exception as e:
        print(f"❌ Неожиданная ошибка: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()