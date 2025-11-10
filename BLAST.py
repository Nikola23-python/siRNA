import pandas as pd
import subprocess
import os
from tqdm import tqdm
from Bio.Seq import Seq
import tempfile

from prepare_rna import df_sense, df_antisense

# Путь к локальной BLAST базе
BLAST_DB = os.path.expanduser("~/blast_dbs/human_refseq_complete")


def check_all_sirna_blast_parameter_local(df_sense, df_antisense, gene_name="ATXN1"):
    """
    ПРОВЕРКА BLAST ПАРАМЕТРА ДЛЯ ВСЕХ siRNA С ЛОКАЛЬНОЙ БАЗОЙ
    """
    print(f"🎯 ПРОВЕРКА BLAST ПАРАМЕТРА ДЛЯ ВСЕХ {len(df_sense)} siRNA")
    print(f"💾 Используется локальная база: {BLAST_DB}")
    print("=" * 70)

    # Проверяем что база существует
    if not os.path.exists(BLAST_DB + ".nhr"):
        print(f"❌ BLAST база не найдена: {BLAST_DB}")
        print("💡 Убедитесь что база создана: ~/blast_dbs/human_refseq_complete")
        return None, 0

    results = []
    total_score = 0

    # Используем tqdm для прогресс-бара
    for idx in tqdm(range(len(df_sense)), desc="Анализ siRNA"):
        sense_row = df_sense.iloc[idx]
        antisense_row = df_antisense.iloc[idx]

        sense_sequence = sense_row['sequence']
        antisense_sequence = antisense_row['sequence']
        fragment_id = sense_row['fragment_id']
        size_nt = sense_row['size_nt']

        # Пропускаем если последовательности не совпадают по длине
        if len(sense_sequence) != len(antisense_sequence):
            continue

        # Проверяем каждую цепь через ЛОКАЛЬНЫЙ BLAST
        sense_score = check_strand_specificity_local(sense_sequence, "СМЫСЛОВАЯ", fragment_id)
        antisense_score = check_strand_specificity_local(antisense_sequence, "АНТИСМЫСЛОВАЯ", fragment_id)

        # Применяем scoring систему из статьи
        sirna_score = calculate_blast_score(sense_score, antisense_score)
        total_score += sirna_score

        # Сохраняем результаты
        results.append({
            'fragment_id': fragment_id,
            'size_nt': size_nt,
            'sense_sequence': sense_sequence,
            'antisense_sequence': antisense_sequence,
            'blast_score': sirna_score,
            'sense_blast_result': sense_score,
            'antisense_blast_result': antisense_score
        })

    # Создаем DataFrame с результатами
    results_df = pd.DataFrame(results)

    # Статистика
    avg_score = total_score / len(results) if results else 0

    print(f"\n📊 ФИНАЛЬНАЯ СТАТИСТИКА:")
    print(f"   Проанализировано siRNA: {len(results)}")
    print(f"   Средний BLAST score: {avg_score:.2f}/2")

    # Анализ распределения scores
    score_distribution = results_df['blast_score'].value_counts().sort_index()
    print(f"   Распределение scores:")
    for score, count in score_distribution.items():
        print(f"     • {score} баллов: {count} siRNA ({count / len(results) * 100:.1f}%)")

    # Лучшие siRNA (score = 2)
    best_sirnas = results_df[results_df['blast_score'] == 2]
    print(f"   Лучших siRNA (2 балла): {len(best_sirnas)}")

    if len(best_sirnas) > 0:
        print(f"   Топ-10 лучших siRNA:")
        for idx, row in best_sirnas.head(10).iterrows():
            print(f"     • {row['fragment_id']} ({row['size_nt']}нт): {row['sense_sequence']}")

    return results_df, avg_score


def check_strand_specificity_local(sequence, strand_type, fragment_id):
    """
    Проверить специфичность одной цепи через ЛОКАЛЬНЫЙ BLAST
    """
    try:
        # Создаем временный файл с последовательностью
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_file:
            temp_file.write(f">{fragment_id}_{strand_type}\n{sequence}\n")
            temp_filename = temp_file.name

        # Выполняем ЛОКАЛЬНЫЙ BLAST с параметрами из статьи
        cmd = [
            "blastn",
            "-query", temp_filename,
            "-db", BLAST_DB,
            "-word_size", "7",  # Статья: word size = 7
            "-evalue", "1000",  # Статья: E-value = 1000-3000
            "-gapopen", "2",  # Статья: gap costs
            "-gapextend", "1",  # Статья: gap costs
            "-outfmt", "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-num_alignments", "10",  # Ограничиваем количество результатов
            "-task", "blastn-short"  # Оптимизировано для коротких последовательностей
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        # Удаляем временный файл
        os.unlink(temp_filename)

        # Анализируем результаты BLAST
        is_good = analyze_blast_results_local(result.stdout, sequence, strand_type, fragment_id)

        return 1 if is_good else 0

    except Exception as e:
        print(f"❌ Ошибка локального BLAST для {fragment_id} ({strand_type}): {e}")
        return 0


def analyze_blast_results_local(blast_output, sequence, strand_type, fragment_id):
    """
    Анализ результатов ЛОКАЛЬНОГО BLAST по критериям из статьи
    """
    lines = [line.strip() for line in blast_output.strip().split('\n') if line.strip()]

    # Если нет выравниваний - отлично!
    if not lines:
        return True

    good_blast = True
    issues = []

    for line in lines:
        parts = line.split(',')
        if len(parts) >= 12:
            try:
                pident = float(parts[2])  # % идентичности
                length = int(parts[3])  # длина выравнивания
                evalue = float(parts[10])  # e-value
                subject_title = parts[12] if len(parts) > 12 else ""  # описание субъекта

                # КРИТЕРИИ из статьи:
                # 1. Покрытие < 78%
                query_coverage = (length / len(sequence)) * 100
                if query_coverage > 78:
                    good_blast = False
                    issues.append(f"покрытие {query_coverage:.1f}%")

                # 2. Совпадений < 15 из 19 (или пропорционально длине)
                matches = int(length * pident / 100)
                max_allowed_matches = min(15, len(sequence) - 2)
                if matches >= max_allowed_matches:
                    good_blast = False
                    issues.append(f"совпадений {matches}/{len(sequence)}")

                # 3. Проверка seed региона (только для смысловой цепи)
                if strand_type == "СМЫСЛОВАЯ" and len(sequence) >= 8:
                    if check_seed_region_issue_local(sequence, subject_title):
                        good_blast = False
                        issues.append("seed регион")

            except (ValueError, IndexError):
                continue

    # Логируем проблемы только для первых нескольких siRNA чтобы не засорять вывод
    if not good_blast and len(issues) > 0 and int(fragment_id.split('_')[1]) < 10:
        print(f"⚠️  {fragment_id} ({strand_type}): {', '.join(set(issues))}")

    return good_blast


def check_seed_region_issue_local(sequence, subject_title):
    """
    Проверка seed региона (позиции 2-8) на off-target эффекты
    """
    if len(sequence) < 8:
        return False

    seed_region = sequence[1:8]  # Позиции 2-8

    # Эвристика: если в описании субъекта нет ATXN1, а seed регион консервативен
    if "ATXN1" not in subject_title.upper():
        # Проверяем консервативность seed региона (высокий GC content)
        gc_count = seed_region.count('G') + seed_region.count('C')
        gc_content = gc_count / len(seed_region)

        # Если seed регион высококонсервативен - возможен off-target
        return gc_content > 0.6

    return False


def calculate_blast_score(sense_score, antisense_score):
    """
    Расчет баллов согласно статье
    """
    if sense_score == 1 and antisense_score == 1:
        return 2
    elif sense_score == 1 or antisense_score == 1:
        return 1
    else:
        return 0


def save_detailed_results(results_df, filename='sirna_blast_detailed_results.csv'):
    """
    Сохранение детальных результатов
    """
    # Добавляем дополнительную информацию
    results_df['has_seed_region'] = results_df['sense_sequence'].apply(
        lambda x: len(x) >= 8
    )

    results_df['seed_sequence'] = results_df['sense_sequence'].apply(
        lambda x: x[1:8] if len(x) >= 8 else 'N/A'
    )

    # Сохраняем
    results_df.to_csv(filename, index=False)
    print(f"💾 Детальные результаты сохранены в '{filename}'")

    return results_df


# 🚀 ЗАПУСК БЫСТРОГО АНАЛИЗА С ЛОКАЛЬНОЙ БАЗОЙ
if __name__ == "__main__":
    print("🚀 ЗАПУСК БЫСТРОГО АНАЛИЗА 1536 siRNA С ЛОКАЛЬНОЙ БАЗОЙ")
    print("⏰ Внимание: теперь это займет МИНУТЫ вместо часов!")
    print("=" * 70)

    try:
        # Запускаем БЫСТРЫЙ анализ с локальной базой
        results_df, avg_score = check_all_sirna_blast_parameter_local(df_sense, df_antisense)

        # Сохраняем детальные результаты
        results_df = save_detailed_results(results_df)

        print(f"\n🎉 АНАЛИЗ ЗАВЕРШЕН ЗА СЕКУНДЫ!")
        print(f"   Итоговый средний score: {avg_score:.2f}/2")
        print(f"   Файл с результатами: 'sirna_blast_detailed_results.csv'")

        # Дополнительная статистика
        print(f"\n📈 ДОПОЛНИТЕЛЬНАЯ СТАТИСТИКА:")
        print(f"   Всего проанализировано: {len(results_df)} siRNA")
        print(f"   siRNA с идеальным score (2): {len(results_df[results_df['blast_score'] == 2])}")
        print(f"   siRNA с хорошим score (1): {len(results_df[results_df['blast_score'] == 1])}")
        print(f"   siRNA с плохим score (0): {len(results_df[results_df['blast_score'] == 0])}")

    except Exception as e:
        print(f"\n❌ Ошибка: {e}")