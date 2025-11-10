import RNA
import pandas as pd

from prepare_rna import df_sense, df_antisense

def sirna_duplex_analysis(df_sense, df_antisense, name=""):
    """Анализ дуплексов для siRNA из DataFrame"""

    results = []

    # Предполагаем, что DataFrame имеют колонку 'sequence'
    for idx, (sense_row, antisense_row) in enumerate(zip(df_sense.iterrows(), df_antisense.iterrows())):
        sense_idx, sense_data = sense_row
        antisense_idx, antisense_data = antisense_row

        sense_seq = sense_data['sequence']  # Получаем последовательность
        antisense_seq = antisense_data['sequence']  # Получаем последовательность

        try:
            # RNA duplex для двух цепей siRNA
            duplex_result = RNA.duplexfold(sense_seq, antisense_seq)

            results.append({
                'sense_id': sense_data.get('fragment_id', f'sense_{idx}'),
                'antisense_id': antisense_data.get('fragment_id', f'antisense_{idx}'),
                'sense_sequence': sense_seq,
                'antisense_sequence': antisense_seq,
                'duplex_structure': duplex_result.structure,
                'duplex_energy': duplex_result.energy,
                'position_i': duplex_result.i,
                'position_j': duplex_result.j
            })

            print(f"Duplex {idx+1}: {duplex_result.energy:.2f} kcal/mol")

        except Exception as e:
            print(f"Ошибка для пары {idx}: {e}")
            results.append({
                'sense_id': sense_data.get('fragment_id', f'sense_{idx}'),
                'antisense_id': antisense_data.get('fragment_id', f'antisense_{idx}'),
                'sense_sequence': sense_seq,
                'antisense_sequence': antisense_seq,
                'duplex_structure': 'ERROR',
                'duplex_energy': None,
                'position_i': None,
                'position_j': None
            })

    return pd.DataFrame(results)

# Запускаем анализ
print("Запускаем RNA duplex...")
duplex_df = sirna_duplex_analysis(df_sense, df_antisense)

# Показываем результаты
print("\nРезультаты RNA duplex:")
print(duplex_df[['sense_id', 'antisense_id', 'duplex_energy', 'duplex_structure']].head())

# Сохраняем в файл
duplex_df.to_csv('rna_duplex_results.csv', index=False)
print(f"\nСохранено результатов: {len(duplex_df)}")
print("Файл: rna_duplex.py_results.csv")