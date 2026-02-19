import pandas as pd
import time
from prepare_rna import editor_rna, ATXN1


def generate_all_fragments_directly():
    """Генерация всех фрагментов напрямую из последовательности"""

    # Параметры
    EXON_START = 16326394
    EXON_END = 16328470
    RNA_START = 781
    RNA_END = 2858
    MIN_LENGTH = 15
    MAX_LENGTH = 30
    MIN_MAF = 0.01

    print(f"{'=' * 80}")
    print("ГЕНЕРАЦИЯ ВСЕХ ФРАГМЕНТОВ 15-30 НТ НАПРЯМУЮ ИЗ ПОСЛЕДОВАТЕЛЬНОСТИ")
    print(f"{'=' * 80}")

    start_time = time.time()

    try:

        # Конвертируем T→U для РНК
        rna_sequence = ATXN1.replace('T', 'U')
        print(f"\n   Длина РНК последовательности: {len(rna_sequence)} нт")

        # Вырезаем интересующий экзон [781:2858]
        exon_sequence = rna_sequence[RNA_START:RNA_END]

        # Проверка длины
        expected_length = EXON_END - EXON_START + 1
        actual_length = len(exon_sequence)


        if actual_length != expected_length:
            print(f"   ⚠ ВНИМАНИЕ: Длины не совпадают! Разница: {abs(actual_length - expected_length)} нт")
            if actual_length < expected_length:
                print(f"   ⚠ Последовательность короче ожидаемой! Возможно неверные координаты.")
        else:
            print(f"   ✓ Длины совпадают!")

        # 2. Загружаем и анализируем SNP
        print(f"\n2. АНАЛИЗ SNP")
        print(f"   {'-' * 40}")

        try:
            gnomad = pd.read_csv("SNP/gnomAD_v4.1.0_ENSG00000124788_2025_12_18_01_08_51.csv")
            print(f"   Загружено {len(gnomad)} строк из gnomAD")

            # Конвертируем типы данных
            gnomad['Position'] = pd.to_numeric(gnomad['Position'], errors='coerce')
            gnomad['Allele Frequency'] = pd.to_numeric(gnomad['Allele Frequency'], errors='coerce')

            # Все SNP в экзоне
            all_snps_in_exon = gnomad[
                (gnomad['Position'] >= EXON_START) &
                (gnomad['Position'] <= EXON_END)
                ]

            print(f"   Всего SNP в экзоне: {len(all_snps_in_exon)}")

            # Частые SNP (MAF ≥ 1%)
            frequent_snps = all_snps_in_exon[
                all_snps_in_exon['Allele Frequency'] >= MIN_MAF
                ].copy()

            print(f"   Частые SNP (MAF ≥ {MIN_MAF * 100}%): {len(frequent_snps)}")

            if len(frequent_snps) == 0:
                print(f"   ⚠ Нет частых SNP в экзоне")
                snp_dict = {}
            else:
                # Создаем словарь SNP: позиция в экзоне → информация о SNP
                snp_dict = {}
                for _, row in frequent_snps.iterrows():
                    genomic_pos = int(row['Position'])

                    # Преобразуем геномную позицию в позицию в экзонной последовательности
                    pos_in_exon = genomic_pos - EXON_START

                    # Проверяем, что позиция в пределах последовательности
                    if 0 <= pos_in_exon < len(exon_sequence):
                        rsid = row.get('rsIDs', f'pos_{genomic_pos}')
                        if pd.isna(rsid) or rsid == '':
                            rsid = f'pos_{genomic_pos}'

                        maf = float(row['Allele Frequency'])

                        if pos_in_exon not in snp_dict:
                            snp_dict[pos_in_exon] = []

                        snp_dict[pos_in_exon].append({
                            'rsid': rsid,
                            'maf': maf,
                            'genomic_pos': genomic_pos,
                            'allele_frequency': maf
                        })

                print(f"   SNP в пределах последовательности: {len(snp_dict)}")

                # Показываем примеры SNP
                if len(snp_dict) > 0:
                    print(f"\n   Примеры частых SNP:")
                    count = 0
                    for pos_in_exon, snp_list in sorted(snp_dict.items())[:5]:
                        for snp in snp_list[:1]:  # Только первый SNP на этой позиции
                            genomic_pos = snp['genomic_pos']
                            nucleotide = exon_sequence[pos_in_exon] if pos_in_exon < len(exon_sequence) else '?'
                            print(
                                f"     Геном {genomic_pos} (РНК {pos_in_exon}='{nucleotide}'): {snp['rsid']} ({snp['maf']:.1%})")
                            count += 1
                            if count >= 3:
                                break
                        if count >= 3:
                            break

        except Exception as e:
            print(f"   ⚠ Ошибка при загрузке SNP: {e}")
            print(f"   Продолжаем без анализа SNP")
            snp_dict = {}


        results = []
        fragment_counter = {}

        # Длина экзонной последовательности
        seq_length = len(exon_sequence)

        print(f"   Длина последовательности: {seq_length} нт")
        print(f"   Диапазон фрагментов: {MIN_LENGTH}-{MAX_LENGTH} нт")

        # Теоретическое количество фрагментов
        total_possible = 0
        for length in range(MIN_LENGTH, MAX_LENGTH + 1):
            total_possible += seq_length - length + 1


        generation_start = time.time()

        # Для каждого возможного размера
        for length in range(MIN_LENGTH, MAX_LENGTH + 1):
            # Количество фрагментов этого размера
            num_fragments_this_size = seq_length - length + 1

            # Для каждой стартовой позиции
            for start_pos in range(num_fragments_this_size):
                end_pos = start_pos + length

                # Последовательность фрагмента
                fragment_seq = exon_sequence[start_pos:end_pos]

                # Геномные координаты
                genomic_start = EXON_START + start_pos
                genomic_end = EXON_START + end_pos - 1

                # Проверяем, что в пределах экзона
                if genomic_end > EXON_END:
                    continue

                # ID фрагмента
                if length not in fragment_counter:
                    fragment_counter[length] = 1
                else:
                    fragment_counter[length] += 1

                frag_id = f"{length}_{fragment_counter[length]}"

                # Поиск SNP во фрагменте
                overlapping_snps = []
                for snp_pos, snp_list in snp_dict.items():
                    # SNP должен быть внутри фрагмента
                    if start_pos <= snp_pos < end_pos:
                        pos_in_fragment = snp_pos - start_pos
                        for snp_info in snp_list:
                            overlapping_snps.append({
                                'rsid': snp_info['rsid'],
                                'pos_in_fragment': pos_in_fragment,
                                'maf': snp_info['maf'],
                                'genomic_pos': snp_info['genomic_pos']
                            })

                # Определяем критические SNP для siRNA
                critical_snp = False
                if 19 <= length <= 23:  # siRNA длины
                    for snp in overlapping_snps:
                        if 2 <= snp['pos_in_fragment'] <= 8:
                            critical_snp = True
                            break

                # Формируем результат
                result = {
                    'fragment_id': frag_id,
                    'sequence': fragment_seq,
                    'length': length,
                    'start_in_exon': start_pos,  # Позиция в экзонной последовательности (0-based)
                    'end_in_exon': end_pos - 1,  # Позиция в экзонной последовательности
                    'genomic_start': genomic_start,
                    'genomic_end': genomic_end,
                    'has_snp': len(overlapping_snps) > 0,
                    'critical_snp': critical_snp,
                    'num_snps': len(overlapping_snps)
                }

                # Добавляем информацию о SNP, если есть
                if overlapping_snps:
                    # Сортируем по MAF (самые частые первыми)
                    overlapping_snps.sort(key=lambda x: -x['maf'])

                    # Форматируем информацию о SNP
                    snp_info_strings = []
                    for snp in overlapping_snps[:5]:  # Берем до 5 самых частых
                        freq_percent = f"{snp['maf'] * 100:.1f}%"
                        snp_info_strings.append(f"{snp['rsid']}(pos{snp['pos_in_fragment']},{freq_percent})")

                    result['top_snp'] = "; ".join(snp_info_strings)
                    result['top_snp_maf'] = overlapping_snps[0]['maf']

                results.append(result)

            # Прогресс для каждого размера
            if num_fragments_this_size > 0:
                elapsed = time.time() - generation_start
                print(f"     {length} нт: {num_fragments_this_size:,} фрагментов ({elapsed:.1f} сек)")

        generation_time = time.time() - generation_start
        print(f"   Генерация завершена за {generation_time:.1f} сек")

        # 4. Создаем DataFrame

        results_df = pd.DataFrame(results)

        # Сортируем по длине и геномной позиции
        results_df = results_df.sort_values(['length', 'genomic_start']).reset_index(drop=True)

        total_fragments = len(results_df)
        total_time = time.time() - start_time

        print(f"   Создано DataFrame: {total_fragments:,} строк × {len(results_df.columns)} столбцов")
        print(f"   Общее время: {total_time:.1f} сек")

        # 5. Статистика
        # Распределение по длинам
        print(f"\nРАСПРЕДЕЛЕНИЕ ПО ДЛИНАМ:")
        for length in range(MIN_LENGTH, MAX_LENGTH + 1):
            count = len(results_df[results_df['length'] == length])
            expected = seq_length - length + 1
            status = "✓" if count == expected else "⚠"

            snp_count = results_df[(results_df['length'] == length) & (results_df['has_snp'])].shape[0]
            percentage = (snp_count / count * 100) if count > 0 else 0

            print(f"  {length:2d} нт: {count:6d} / {expected:6d} {status} (с SNP: {snp_count:5d}, {percentage:5.1f}%)")

        # Статистика по SNP
        fragments_with_snp = results_df['has_snp'].sum()
        percentage_with_snp = (fragments_with_snp / total_fragments * 100) if total_fragments > 0 else 0

        print(f"\nСТАТИСТИКА ПО SNP:")
        print(f"  Фрагментов с SNP: {fragments_with_snp:,} ({percentage_with_snp:.1f}%)")

        # siRNA статистика
        sirna_fragments = results_df[results_df['length'].between(19, 23)]
        sirna_total = len(sirna_fragments)

        if sirna_total > 0:
            sirna_without_snp = len(sirna_fragments[~sirna_fragments['has_snp']])
            sirna_critical = len(sirna_fragments[sirna_fragments['critical_snp']])

            sirna_no_snp_perc = (sirna_without_snp / sirna_total * 100) if sirna_total > 0 else 0
            sirna_critical_perc = (sirna_critical / sirna_total * 100) if sirna_total > 0 else 0


        # Основной файл
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        main_file = f'SNP/all_fragments_{MIN_LENGTH}_{MAX_LENGTH}_with_snps_{timestamp}.csv'
        results_df.to_csv(main_file, index=False)
        print(f"✓ Все фрагменты сохранены в: {main_file}")
        print(f"  Размер: {total_fragments:,} строк × {len(results_df.columns)} столбцов")



        # Фрагменты с SNP
        if fragments_with_snp > 0:
            snp_file = f'SNP/fragments_with_snps_{timestamp}.csv'
            fragments_with_snps = results_df[results_df['has_snp']].copy()
            fragments_with_snps.to_csv(snp_file, index=False)
            print(f"✓ Фрагменты с SNP сохранены в: {snp_file}")
            print(f"  Количество: {fragments_with_snp:,}")


        return results_df

    except Exception as e:
        print(f"\n{'=' * 80}")
        print("ОШИБКА ВЫПОЛНЕНИЯ")
        print(f"{'=' * 80}")
        print(f"Тип ошибки: {type(e).__name__}")
        print(f"Сообщение: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    print(f"НАЧАЛО АНАЛИЗА: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'=' * 80}")

    results = generate_all_fragments_directly()

    if results is not None:
        print(f"\nЗАВЕРШЕНО: {time.strftime('%Y-%m-%d %H:%M:%S')}")

        if len(results) > 0:
            total = len(results)
            with_snp = results['has_snp'].sum()
            sirna = results[results['length'].between(19, 23)]
            sirna_no_snp = len(sirna[~sirna['has_snp']]) if len(sirna) > 0 else 0

            print(f"• Всего фрагментов: {total:,}")
            print(f"• Фрагментов с SNP: {with_snp:,} ({with_snp / total * 100:.1f}%)")

            # Показать диапазон координат
            print(f"\n• Диапазон геномных координат: {results['genomic_start'].min()} - {results['genomic_end'].max()}")
            print(f"• Ожидаемый диапазон: 16326394 - 16328470")
    else:
        print(f"\nАнализ завершен с ошибкой")