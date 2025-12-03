import pandas as pd
import subprocess
import os
import tempfile
from collections import defaultdict
from tqdm import tqdm
import multiprocessing as mp
from functools import partial

# Пути к данным
from prepare_rna import df_sense, df_antisense

# Путь к BLAST базе
BLAST_DB = "/home/nikolay/blast_dbs/human_refseq_complete"


def collect_all_unique_sequences():
    print("📊 СБОР ВСЕХ УНИКАЛЬНЫХ ПОСЛЕДОВАТЕЛЬНОСТЕЙ")

    unique_sequences = set()
    seq_to_sirna = defaultdict(list)  # Сопоставление последовательности с siRNA

    # Собираем sense strands
    for idx, row in df_sense.iterrows():
        seq = row['sequence']
        fragment_id = row['fragment_id']
        unique_sequences.add(seq)
        seq_to_sirna[seq].append(('sense', fragment_id))

    # Собираем antisense strands
    for idx, row in df_antisense.iterrows():
        seq = row['sequence']
        fragment_id = row['fragment_id']
        unique_sequences.add(seq)
        seq_to_sirna[seq].append(('antisense', fragment_id))

    print(f"   Всего последовательностей: {len(df_sense) + len(df_antisense)}")
    print(f"   Уникальных последовательностей: {len(unique_sequences)}")
    print(f"   Экономия: {100 - len(unique_sequences) / (len(df_sense) + len(df_antisense)) * 100:.1f}%")

    return list(unique_sequences), seq_to_sirna


def create_batch_files(sequences, batch_size=1000):
    """
    Создание батч-файлов для BLAST
    """
    print(f"📁 СОЗДАНИЕ БАТЧ-ФАЙЛОВ (размер батча: {batch_size})")

    batches = []
    batch_dir = "blast_batches"
    os.makedirs(batch_dir, exist_ok=True)

    for i in range(0, len(sequences), batch_size):
        batch_seqs = sequences[i:i + batch_size]
        batch_file = os.path.join(batch_dir, f"batch_{i // batch_size}.fasta")

        with open(batch_file, 'w') as f:
            for idx, seq in enumerate(batch_seqs):
                seq_id = f"seq_{i + idx}"
                f.write(f">{seq_id}\n{seq}\n")

        batches.append({
            'file': batch_file,
            'sequences': batch_seqs,
            'start_idx': i,
            'batch_num': i // batch_size
        })

    print(f"   Создано {len(batches)} батч-файлов")
    return batches


def run_batch_blast(batch_info):
    """
    Запуск BLAST для одного батча
    """
    batch_file = batch_info['file']
    batch_num = batch_info['batch_num']
    batch_size = len(batch_info['sequences'])

    # Оптимизированные параметры для быстрого BLAST
    cmd = [
        "blastn",
        "-query", batch_file,
        "-db", BLAST_DB,
        "-task", "blastn-short",
        "-word_size", "11",  # Увеличили для скорости
        "-evalue", "10",  # Быстрая фильтрация
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        "-num_alignments", "3",  # Только топ-3 результата
        "-max_hsps", "1",  # Только лучшее выравнивание
        "-perc_identity", "80",  # Минимум 80% идентичности
        "-qcov_hsp_perc", "80",  # Минимум 80% покрытия
        "-dust", "yes",  # Фильтр низкокомплексных регионов
        "-num_threads", "2"  # Используем 2 потока
    ]

    try:
        # Запускаем BLAST
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)  # 5 минут таймаут

        # Парсим результаты
        blast_results = defaultdict(list)
        for line in result.stdout.strip().split('\n'):
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    seq_id = parts[0]  # Например: seq_0
                    blast_results[seq_id].append(line)

        # Логируем прогресс
        with open(f"blast_batches/batch_{batch_num}_log.txt", 'w') as log:
            log.write(f"Батч {batch_num}: {batch_size} последовательностей\n")
            log.write(f"Найдено совпадений: {len(blast_results)}\n")
            if result.stderr:
                log.write(f"Ошибки: {result.stderr}\n")

        return {
            'batch_num': batch_num,
            'blast_results': dict(blast_results),
            'status': 'success',
            'total_sequences': batch_size,
            'matches_found': len(blast_results)
        }

    except subprocess.TimeoutExpired:
        return {
            'batch_num': batch_num,
            'blast_results': {},
            'status': 'timeout',
            'total_sequences': batch_size,
            'matches_found': 0
        }
    except Exception as e:
        return {
            'batch_num': batch_num,
            'blast_results': {},
            'status': f'error: {str(e)}',
            'total_sequences': batch_size,
            'matches_found': 0
        }


def process_all_batches_parallel(batches, num_workers=4):
    """
    Параллельная обработка всех батчей
    """
    print(f"⚡ ПАРАЛЛЕЛЬНАЯ ОБРАБОТКА BLAST ({num_workers} потоков)")

    with mp.Pool(processes=num_workers) as pool:
        results = list(tqdm(
            pool.imap(run_batch_blast, batches),
            total=len(batches),
            desc="BLAST обработка"
        ))

    return results


def analyze_blast_results_simple(blast_output, sequence):
    """
    Простой анализ результатов BLAST:
    Возвращает True если специфично (нет off-target), False если есть проблемы
    """
    if not blast_output:
        return True, "No hits"  # Нет совпадений - отлично!

    lines = blast_output.strip().split('\n')

    for line in lines:
        parts = line.split('\t')
        if len(parts) >= 13:
            subject_title = parts[12]
            pident = float(parts[2])
            length = int(parts[3])

            # Если это не наш целевой ген ATXN1
            if "ATXN1" not in subject_title.upper():
                # Проверяем качество совпадения
                coverage = (length / len(sequence)) * 100
                if coverage > 70 and pident > 70:  # Хорошее совпадение с другим геном
                    return False, f"Match to {subject_title} ({coverage:.1f}%, {pident:.1f}% id)"

    return True, "Specific"


def compile_results(batch_results, sequences, seq_to_sirna):
    """
    Компиляция всех результатов
    """
    print("📊 КОМПИЛЯЦИЯ РЕЗУЛЬТАТОВ")

    # Создаем словарь для результатов каждой последовательности
    sequence_results = {}

    # Обрабатываем результаты батчей
    for batch in batch_results:
        if batch['status'] != 'success':
            continue

        for seq_id, blast_output in batch['blast_results'].items():
            # Извлекаем индекс из seq_id (например, seq_123 -> 123)
            idx = int(seq_id.split('_')[1])

            # Находим соответствующую последовательность
            if idx < len(sequences):
                seq = sequences[idx]
                is_specific, reason = analyze_blast_results_simple('\n'.join(blast_output), seq)
                sequence_results[seq] = {
                    'specific': is_specific,
                    'reason': reason,
                    'hits_count': len(blast_output)
                }

    # Теперь собираем результаты по siRNA
    sirna_results = []

    print("🧬 СБОР РЕЗУЛЬТАТОВ ПО siRNA...")
    for fragment_id in tqdm(df_sense['fragment_id'].unique(), desc="Обработка siRNA"):
        # Находим sense и antisense последовательности для этой siRNA
        sense_row = df_sense[df_sense['fragment_id'] == fragment_id].iloc[0]
        anti_row = df_antisense[df_antisense['fragment_id'] == fragment_id].iloc[0]

        sense_seq = sense_row['sequence']
        anti_seq = anti_row['sequence']
        size_nt = sense_row['size_nt']

        # Получаем результаты для sense
        sense_result = sequence_results.get(sense_seq, {'specific': True, 'reason': 'No data', 'hits_count': 0})
        anti_result = sequence_results.get(anti_seq, {'specific': True, 'reason': 'No data', 'hits_count': 0})

        # Подсчет BLAST score
        if sense_result['specific'] and anti_result['specific']:
            blast_score = 2
        elif sense_result['specific'] or anti_result['specific']:
            blast_score = 1
        else:
            blast_score = 0

        sirna_results.append({
            'fragment_id': fragment_id,
            'size_nt': size_nt,
            'sense_sequence': sense_seq,
            'antisense_sequence': anti_seq,
            'sense_specific': sense_result['specific'],
            'antisense_specific': anti_result['specific'],
            'sense_hits': sense_result['hits_count'],
            'antisense_hits': anti_result['hits_count'],
            'sense_reason': sense_result['reason'],
            'antisense_reason': anti_result['reason'],
            'blast_score': blast_score
        })

    return pd.DataFrame(sirna_results)


def save_and_analyze_results(results_df):
    """
    Сохранение и анализ результатов
    """
    print("\n💾 СОХРАНЕНИЕ РЕЗУЛЬТАТОВ")

    # Сохраняем все результаты
    results_df.to_csv('sirna_blast_only_results.csv', index=False)

    # Фильтруем только хорошие siRNA (score 2)
    good_sirnas = results_df[results_df['blast_score'] == 2]
    good_sirnas.to_csv('sirna_blast_good_results.csv', index=False)

    # Статистика
    total = len(results_df)
    score_2 = len(good_sirnas)
    score_1 = len(results_df[results_df['blast_score'] == 1])
    score_0 = len(results_df[results_df['blast_score'] == 0])

    print(f"📊 СТАТИСТИКА BLAST ПРОВЕРКИ:")
    print(f"   Всего siRNA: {total}")
    print(f"   Score 2 (обе цепи специфичны): {score_2} ({score_2 / total * 100:.1f}%)")
    print(f"   Score 1 (одна цепь специфична): {score_1} ({score_1 / total * 100:.1f}%)")
    print(f"   Score 0 (обе цепи неспецифичны): {score_0} ({score_0 / total * 100:.1f}%)")

    # Топ-10 лучших siRNA
    print(f"\n🏆 ТОП-10 ЛУЧШИХ siRNA (по BLAST score):")
    for idx, row in good_sirnas.head(10).iterrows():
        print(f"   {row['fragment_id']} ({row['size_nt']}нт):")
        print(f"     Sense: {row['sense_sequence']}")
        print(f"     Anti:  {row['antisense_sequence']}")

    # Причины проблем
    if score_0 > 0:
        problematic = results_df[results_df['blast_score'] == 0]
        print(f"\n⚠️  ПРИЧИНЫ ПРОБЛЕМ (первые 5):")
        for idx, row in problematic.head(5).iterrows():
            print(f"   {row['fragment_id']}:")
            if not row['sense_specific']:
                print(f"     Sense: {row['sense_reason']}")
            if not row['antisense_specific']:
                print(f"     Anti:  {row['antisense_reason']}")

    print(f"\n💾 ФАЙЛЫ:")
    print(f"   • Все результаты: sirna_blast_only_results.csv")
    print(f"   • Хорошие siRNA (score 2): sirna_blast_good_results.csv")


def main_full_blast_check():
    """
    Полная проверка ВСЕХ siRNA через BLAST
    """
    print("=" * 80)
    print("🔥 ПОЛНАЯ BLAST ПРОВЕРКА ВСЕХ siRNA (32888 пар)")
    print("=" * 80)

    # Шаг 1: Сбор всех уникальных последовательностей
    sequences, seq_to_sirna = collect_all_unique_sequences()

    # Шаг 2: Создание батч-файлов
    batches = create_batch_files(sequences, batch_size=500)  # 500 последовательностей в батче

    # Шаг 3: Параллельная обработка
    print("\n⚡ ЗАПУСК ПАРАЛЛЕЛЬНОГО BLAST...")
    print("   Это может занять 30-60 минут в зависимости от системы")

    batch_results = process_all_batches_parallel(batches, num_workers=4)

    # Шаг 4: Анализ статусов батчей
    print("\n📈 СТАТУСЫ БАТЧЕЙ:")
    status_counts = {}
    for batch in batch_results:
        status = batch['status']
        status_counts[status] = status_counts.get(status, 0) + 1

    for status, count in status_counts.items():
        print(f"   {status}: {count} батчей")

    # Шаг 5: Компиляция результатов
    print("\n📊 КОМПИЛЯЦИЯ ВСЕХ РЕЗУЛЬТАТОВ...")
    results_df = compile_results(batch_results, sequences, seq_to_sirna)

    # Шаг 6: Сохранение и анализ
    save_and_analyze_results(results_df)

    print("\n✅ ПОЛНАЯ BLAST ПРОВЕРКА ЗАВЕРШЕНА!")
    print(f"   Проверено: {len(results_df)} siRNA")
    print(f"   Найдено специфичных: {len(results_df[results_df['blast_score'] == 2])}")

    return results_df


def quick_blast_check():
    """
    Быстрая проверка (только первые 1000 siRNA)
    """
    print("🚀 БЫСТРАЯ ПРОВЕРКА (первые 1000 siRNA)")

    # Берем первые 1000 siRNA
    df_sense_small = df_sense.head(1000)
    df_antisense_small = df_antisense.head(1000)

    # Собираем уникальные последовательности
    sequences = set()
    for seq in df_sense_small['sequence']:
        sequences.add(seq)
    for seq in df_antisense_small['sequence']:
        sequences.add(seq)

    sequences = list(sequences)
    print(f"   Уникальных последовательностей: {len(sequences)}")

    # Создаем один батч
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        for idx, seq in enumerate(sequences):
            f.write(f">seq_{idx}\n{seq}\n")
        batch_file = f.name

    # Запускаем BLAST
    cmd = [
        "blastn",
        "-query", batch_file,
        "-db", BLAST_DB,
        "-task", "blastn-short",
        "-word_size", "11",
        "-evalue", "10",
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
        "-num_alignments", "3",
        "-max_hsps", "1",
        "-perc_identity", "70",
        "-qcov_hsp_perc", "50",
        "-dust", "yes"
    ]

    print("   Запуск BLAST...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    os.unlink(batch_file)

    # Парсим результаты
    blast_results = defaultdict(list)
    for line in result.stdout.strip().split('\n'):
        if line and not line.startswith('#'):
            parts = line.split('\t')
            if len(parts) >= 2:
                seq_id = parts[0]
                blast_results[seq_id].append(line)

    print(f"   Найдено совпадений для {len(blast_results)} последовательностей")

    # Анализируем результаты
    results = []
    for idx in range(len(df_sense_small)):
        sense_row = df_sense_small.iloc[idx]
        anti_row = df_antisense_small.iloc[idx]

        sense_seq = sense_row['sequence']
        anti_seq = anti_row['sequence']

        # Находим индексы последовательностей
        sense_idx = sequences.index(sense_seq) if sense_seq in sequences else -1
        anti_idx = sequences.index(anti_seq) if anti_seq in sequences else -1

        # Проверяем специфичность
        sense_specific = True
        anti_specific = True

        if sense_idx != -1 and f"seq_{sense_idx}" in blast_results:
            sense_output = '\n'.join(blast_results[f"seq_{sense_idx}"])
            sense_specific, _ = analyze_blast_results_simple(sense_output, sense_seq)

        if anti_idx != -1 and f"seq_{anti_idx}" in blast_results:
            anti_output = '\n'.join(blast_results[f"seq_{anti_idx}"])
            anti_specific, _ = analyze_blast_results_simple(anti_output, anti_seq)

        # Подсчет score
        if sense_specific and anti_specific:
            blast_score = 2
        elif sense_specific or anti_specific:
            blast_score = 1
        else:
            blast_score = 0

        results.append({
            'fragment_id': sense_row['fragment_id'],
            'size_nt': sense_row['size_nt'],
            'sense_sequence': sense_seq,
            'antisense_sequence': anti_seq,
            'blast_score': blast_score
        })

    results_df = pd.DataFrame(results)

    # Статистика
    score_2 = len(results_df[results_df['blast_score'] == 2])
    score_1 = len(results_df[results_df['blast_score'] == 1])
    score_0 = len(results_df[results_df['blast_score'] == 0])

    print(f"\n📊 СТАТИСТИКА:")
    print(f"   Score 2: {score_2} siRNA")
    print(f"   Score 1: {score_1} siRNA")
    print(f"   Score 0: {score_0} siRNA")

    # Сохраняем
    results_df.to_csv('sirna_blast_quick_check.csv', index=False)
    print(f"\n💾 Результаты сохранены в 'sirna_blast_quick_check.csv'")

    return results_df


if __name__ == "__main__":
    print("Выберите режим:")
    print("1. Полная проверка всех siRNA (32888 пар) - 30-60 минут")
    print("2. Быстрая проверка (первые 1000 siRNA) - 5 минут")

    choice = input("Введите 1 или 2: ")

    if choice == "1":
        results = main_full_blast_check()
    elif choice == "2":
        results = quick_blast_check()
    else:
        print("Неверный выбор. Запускаю быструю проверку...")
        results = quick_blast_check()