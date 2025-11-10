import pandas as pd
import RNA
from prepare_rna import df_sense, df_antisense

# Применяем RNA fold ко всем последовательностям
def apply_rna_fold(df_sense):
    """Применяет RNA fold к DataFrame и возвращает результаты"""

    results = []

    for idx, row in df_sense.iterrows():
        fragment_id = row['fragment_id']
        sequence = row['sequence']

        try:
            # ПРИМЕНЯЕМ RNA FOLD - основная функция
            structure, mfe = RNA.fold(sequence)

            results.append({
                'fragment_id': fragment_id,
                'sequence': sequence,
                'structure': structure,  # Складчатая структура в скобочной нотации
                'mfe': mfe,  # Minimal Free Energy
                'sequence_length': len(sequence)
            })

        except Exception as e:
            print(f"Ошибка для {fragment_id}: {e}")
            results.append({
                'fragment_id': fragment_id,
                'sequence': sequence,
                'structure': 'ERROR',
                'mfe': None,
                'sequence_length': len(sequence)
            })

    return pd.DataFrame(results)


# Запускаем RNA fold
print("Запускаем RNA fold...")
folded_df_sense = apply_rna_fold(df_sense)

# Показываем результаты
print("\nРезультаты RNA fold:")
print(folded_df_sense[['fragment_id', 'sequence', 'structure', 'mfe']].head(10))
# Проверяем различные типы структур
print("Примеры структур:")
for i, row in folded_df_sense.head(5).iterrows():
    print(f"{row['fragment_id']}: {row['structure']} (энергия: {row['mfe']:.2f})")

# Статистика по энергиям
print(f"\nСтатистика по свободной энергии:")
print(f"Минимальная: {folded_df_sense['mfe'].min():.2f}")
print(f"Максимальная: {folded_df_sense['mfe'].max():.2f}")
print(f"Средняя: {folded_df_sense['mfe'].mean():.2f}")

# Сохраняем в файл
folded_df_sense.to_csv('rna_fold_results.csv', index=False)
print(f"\nСохранено результатов: {len(folded_df_sense)}")
print("Файл: rna_fold_results.csv")