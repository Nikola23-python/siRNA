import time
from Bio.Blast import NCBIWWW, NCBIXML
from prepare_rna import df_sense, df_antisense
import pandas as pd
from typing import Dict, Any

import time
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
from typing import Dict, Any

# Импортируем только если выполняется как основной скрипт
try:
    from prepare_rna import df_sense, df_antisense
except ImportError:
    # Для тестирования создадим пустые DataFrame
    df_sense = pd.DataFrame()
    df_antisense = pd.DataFrame()


class BlastScorer:
    def __init__(self, max_coverage: float = 0.78, max_homology: int = 7):
        """
        Инициализация скорринговой системы

        Args:
            max_coverage (float): максимально допустимое покрытие (0.78 = 78%)
            max_homology (int): максимально допустимая гомология в нуклеотидах
        """
        self.max_coverage = max_coverage
        self.max_homology = max_homology

    def run_blast_analysis(self, sequences_df, database="nt", wait_time=5):
        """
        Выполняет BLAST анализ для последовательностей в DataFrame
        """
        results = []

        for idx, row in sequences_df.iterrows():
            fragment_id = row['fragment_id']
            sequence = row['sequence']

            print(f"Анализируем {fragment_id}: {sequence}")

            try:
                # Выполняем BLAST поиск
                result_handle = NCBIWWW.qblast("blastn", database, sequence)

                # Парсим результаты
                blast_records = NCBIXML.parse(result_handle)

                max_coverage = 0
                max_homology = 0

                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            # Рассчитываем покрытие
                            coverage = hsp.align_length / len(sequence)
                            max_coverage = max(max_coverage, coverage)

                            # Ищем максимальную гомологию
                            max_homology = max(max_homology, hsp.align_length)

                results.append({
                    'fragment_id': fragment_id,
                    'max_coverage': max_coverage,
                    'max_homology': max_homology
                })

                # Закрываем handle для освобождения ресурсов
                result_handle.close()

                # Пауза чтобы не перегружать сервер NCBI
                time.sleep(wait_time)

            except Exception as e:
                print(f"Ошибка для {fragment_id}: {e}")
                results.append({
                    'fragment_id': fragment_id,
                    'max_coverage': 1.0,  # В случае ошибки считаем худший вариант
                    'max_homology': len(sequence)
                })

        return pd.DataFrame(results)

    def evaluate_blast_dataframes(self, df_sense: pd.DataFrame, df_antisense: pd.DataFrame,
                                  id_column: str = "fragment_id") -> pd.DataFrame:
        """
        Оценивает BLAST результаты из двух DataFrame и возвращает один DataFrame с баллами

        Args:
            df_sense (pd.DataFrame): DataFrame с BLAST результатами для смысловых цепей
            df_antisense (pd.DataFrame): DataFrame с BLAST результатами для антисмысловых цепей
            id_column (str): название колонки с идентификатором последовательности

        Returns:
            pd.DataFrame: DataFrame с результатами оценки
        """
        # Проверяем, что ID колонка существует в обоих DataFrame
        if id_column not in df_sense.columns or id_column not in df_antisense.columns:
            raise ValueError(f"Колонка {id_column} не найдена в одном из DataFrame")

        # Объединяем данные по ID
        merged_df = pd.merge(df_sense, df_antisense, on=id_column,
                             suffixes=('_sense', '_antisense'))

        # Применяем функцию оценки к каждой строке
        results = []
        for _, row in merged_df.iterrows():
            result = self._evaluate_row(row, id_column)
            results.append(result)

        # Создаем итоговый DataFrame
        results_df = pd.DataFrame(results)
        return results_df

    def _evaluate_row(self, row: pd.Series, id_column: str) -> Dict[str, Any]:
        """
        Оценивает одну строку с данными смысловой и антисмысловой цепи
        """
        # Подготавливаем данные для смысловой цепи
        sense_data = {
            "max_coverage": row.get("max_coverage_sense", row.get("coverage_sense", 0)),
            "max_homology": row.get("max_homology_sense", row.get("homology_sense", 0))
        }

        # Подготавливаем данные для антисмысловой цепи
        antisense_data = {
            "max_coverage": row.get("max_coverage_antisense", row.get("coverage_antisense", 0)),
            "max_homology": row.get("max_homology_antisense", row.get("homology_antisense", 0))
        }

        # Оцениваем каждую цепь отдельно
        sense_score = self._evaluate_single_strand(sense_data)
        antisense_score = self._evaluate_single_strand(antisense_data)

        # Рассчитываем общий балл по правилам
        total_score = 0
        if sense_score["passed"] and antisense_score["passed"]:
            total_score = 2
        elif sense_score["passed"] or antisense_score["passed"]:
            total_score = 1

        # Формируем результат
        return {
            id_column: row[id_column],
            "blast_total_score": total_score,
            "sense_passed": sense_score["passed"],
            "sense_coverage": sense_score["coverage_value"],
            "sense_homology": sense_score["max_homology_found"],
            "sense_coverage_check": sense_score["coverage_check"],
            "sense_homology_check": sense_score["homology_check"],
            "antisense_passed": antisense_score["passed"],
            "antisense_coverage": antisense_score["coverage_value"],
            "antisense_homology": antisense_score["max_homology_found"],
            "antisense_coverage_check": antisense_score["coverage_check"],
            "antisense_homology_check": antisense_score["homology_check"],
            "explanation": self._generate_explanation(sense_score, antisense_score, total_score)
        }

    def _evaluate_single_strand(self, blast_data: Dict) -> Dict[str, Any]:
        """
        Оценивает одну цепь siRNA
        """
        # Проверяем покрытие
        coverage_ok = self._check_coverage(blast_data)

        # Проверяем гомологию
        homology_ok = self._check_homology(blast_data)

        # Цепь проходит, если обе проверки пройдены
        passed = coverage_ok and homology_ok

        return {
            "passed": passed,
            "coverage_check": coverage_ok,
            "homology_check": homology_ok,
            "coverage_value": blast_data.get("max_coverage", 0),
            "max_homology_found": blast_data.get("max_homology", 0)
        }

    def _check_coverage(self, blast_data: Dict) -> bool:
        """Проверяет покрытие query coverage"""
        max_coverage = blast_data.get("max_coverage", 0)
        return max_coverage < self.max_coverage

    def _check_homology(self, blast_data: Dict) -> bool:
        """Проверяет гомологию по длине последовательности"""
        max_homology = blast_data.get("max_homology", 0)
        return max_homology < self.max_homology

    def _generate_explanation(self, sense_score: Dict, antisense_score: Dict, total_score: int) -> str:
        """Генерирует текстовое объяснение результата"""
        explanations = []

        if total_score == 2:
            explanations.append("ОБЕ цепи показали отличные результаты BLAST")
        elif total_score == 1:
            passed_strands = []
            if sense_score["passed"]:
                passed_strands.append("смысловая")
            if antisense_score["passed"]:
                passed_strands.append("антисмысловая")
            explanations.append(f"Только {', '.join(passed_strands)} цепь прошла проверку")
        else:
            explanations.append("Ни одна цепь не прошла проверку BLAST")

        # Добавляем детали по каждой цепи
        for strand_type, score in [("Смысловая", sense_score), ("Антисмысловая", antisense_score)]:
            if not score["passed"]:
                reasons = []
                if not score["coverage_check"]:
                    reasons.append(f"покрытие {score['coverage_value']:.1%} > {self.max_coverage:.1%}")
                if not score["homology_check"]:
                    reasons.append(f"гомология {score['max_homology_found']} нт > {self.max_homology} нт")
                if reasons:
                    explanations.append(f"{strand_type} цепь: {', '.join(reasons)}")

        return "; ".join(explanations)

    def get_scoring_summary(self, results_df: pd.DataFrame) -> pd.DataFrame:
        """
        Возвращает сводку по результатам scoring'а

        Args:
            results_df (pd.DataFrame): DataFrame с результатами оценки

        Returns:
            pd.DataFrame: Сводная статистика
        """
        if len(results_df) == 0:
            return pd.DataFrame([{
                'total_sequences': 0,
                'score_2_count': 0,
                'score_1_count': 0,
                'score_0_count': 0,
                'score_2_percent': 0,
                'sense_pass_rate': 0,
                'antisense_pass_rate': 0
            }])

        summary = {
            'total_sequences': len(results_df),
            'score_2_count': len(results_df[results_df['blast_total_score'] == 2]),
            'score_1_count': len(results_df[results_df['blast_total_score'] == 1]),
            'score_0_count': len(results_df[results_df['blast_total_score'] == 0]),
            'score_2_percent': len(results_df[results_df['blast_total_score'] == 2]) / len(results_df) * 100,
            'sense_pass_rate': results_df['sense_passed'].mean() * 100,
            'antisense_pass_rate': results_df['antisense_passed'].mean() * 100
        }

        return pd.DataFrame([summary])


def main():
    """Основная функция для запуска BLAST анализа"""

    # Проверяем, что данные загружены
    if df_sense.empty or df_antisense.empty:
        print("Ошибка: DataFrame с последовательностями пусты!")
        return

    # Создаем оценщик
    scorer = BlastScorer(max_coverage=0.78, max_homology=7)

    # Шаг 2: Выполняем BLAST анализ (это может занять время!)
    print("Запускаем BLAST анализ для смысловых цепей...")
    blast_sense = scorer.run_blast_analysis(df_sense)
    print("Запускаем BLAST анализ для антисмысловых цепей...")
    blast_antisense = scorer.run_blast_analysis(df_antisense)

    # Шаг 3: Объединяем BLAST результаты с исходными данными
    df_sense_with_blast = pd.merge(df_sense, blast_sense, on='fragment_id')
    df_antisense_with_blast = pd.merge(df_antisense, blast_antisense, on='fragment_id')

    # Запускаем оценку
    results_df = scorer.evaluate_blast_dataframes(df_sense_with_blast, df_antisense_with_blast)

    print("\nРезультаты оценки:")
    print(results_df[['fragment_id', 'blast_total_score', 'sense_passed', 'antisense_passed', 'explanation']])

    # Получаем сводку
    summary = scorer.get_scoring_summary(results_df)
    print("\nСводка по результатам:")
    print(summary)
    summary.to_csv('blast_scoring_summary.csv', index=False, encoding='utf-8')
    print("✓ Сводка сохранена в blast_scoring_summary.csv")
    return results_df


if __name__ == "__main__":
    results = main()