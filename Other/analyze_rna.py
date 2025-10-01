# from datetime import datetime
# import pandas as pd
# from Bio import SeqIO
# sequence1 = SeqIO.read("sequence.fasta", "fasta")
# ATXN1 = str(sequence1.seq)
#
#
# class editor_rna:
#     """Метод класса для работы с размерами RNA"""
#     def __init__(self):
#         self.sequence = ATXN1
#         self.length = len(self.sequence)
#         self._fragments = []
#         self._gc_content = []
#
#     def info_rna(self):
#         """Метод выводит последовательность РНК"""
#         print(self.sequence)
#
#     def slice_rna(self, start_size=15, end_size=30, filename=None):
#         # start = int(input("Введите начальный индекс: "))
#         # end = int(input("Введите конечный индекс: "))
#         # frag = self.sequence[start:end]
#         frag = self.sequence[781:2858]
#         # start_size = int(input("Введите начальный размер: "))
#         # end_size = int(input("Введите конечный размер: ")) + 1
#         start_size = 15
#         end_size = 30 + 1
#
#         self._fragments = []
#
#         print("🔪 Нарезка РНК последовательности...")
#
#         for fragment_size in range(start_size, end_size):
#             for start_pos in range(0, len(frag), fragment_size):
#                 end_pos = start_pos + fragment_size
#                 if end_pos > len(frag):
#                     break
#                 fragment = frag[start_pos:start_pos + fragment_size]
#                 self._fragments.append(
#                     # 'fragment_id': f"{fragment_size}_{start_pos + 1}",
#                     # 'size_nt': fragment_size,
#                     fragment
#                     # 'sequence_length': len(frag)
#                 )
#         fragments = str(self._fragments)
#         return fragments
#
# b = editor_rna()
# print(b.slice_rna())
#
# class analyz_rna:
#     """Конкретные методы класса для обработки RNA"""
#
#     def __init__(self):
#         self.parts = None
#         self._gc_content = None
#
#     #     self.editor_rna = editor_rna()  # 🔹 Создаем другой класс
#     #
#     #     # 🔹 ПЕРЕМЕННЫЕ ДЛЯ ХРАНЕНИЯ ДАННЫХ ИЗ ДРУГОГО КЛАССА
#     #     self.parts = []  # Сюда скопируем fragments_list
#     #
#     # def import_fragments(self):
#     #     """Копируем переменную из другого класса"""
#     #     # Генерируем фрагменты в другом классе
#     #     self.editor_rna.slice_rna(self.sequence)
#     #     # 🔹 КОПИРУЕМ ПЕРЕМЕННУЮ ИЗ ДРУГОГО КЛАССА
#     #     self.parts = self.editor_rna._fragments# 🔹 Прямое копирование
#     #
#     #     print(f"✅ Импортировано фрагментов: {len(self.parts)}")
#     #     return self.parts
#
#     def analyze_gc(self, start_size=15, end_size=30, filename=None):
#         """Анализирует GC-content для списка последовательностей и возвращает DataFrame"""
#         self.parts = editor_rna.slice_rna()
#         self._gc_content = []  # Инициализируем пустой список
#         for i, item in enumerate(self.parts):
#             # Определяем, что именно является последовательностью
#             if isinstance(item, tuple):
#                 # Если это кортеж, берем первый элемент (предположительно последовательность)
#                 fragment = item[0]
#             elif isinstance(item, str):
#                 # Если это строка, используем как есть
#                 fragment = item
#             else:
#                 # Для других типов преобразуем в строку
#                 fragment = str(item)
#
#             # Проверяем, что fragment - строка
#             if not isinstance(fragment, str):
#                 fragment = str(fragment)
#             length = len(fragment)
#             gc_count = fragment.upper().count('G') + fragment.upper().count('C')
#             total_length = length
#             if total_length > 0:
#                 gc_percent = (gc_count / total_length) * 100
#             else:
#                 gc_percent = 0.0
#             at_count = length - gc_count
#
#             # ДОБАВЛЯЕМ в список, а не перезаписываем
#             self._gc_content.append({
#                 'sequence_id': f'seq_{1:03d}',
#                 'sequence_preview': fragment,
#                 'length': length,
#                 'gc_count': gc_count,
#                 'at_count': at_count,
#                 'gc_percent': gc_percent,
#                 'gc_category': 'Low' if gc_percent < 36 else 'High' if gc_percent > 52 else 'Optimal'
#             })
#         return pd.DataFrame(self._gc_content)
#
#             # df = pd.DataFrame(results)
#             # if filename is None:
#             #     timestamp = datetime.now().strftime("%Y%m%d_%H%M")
#             #     filename = f"rna_fragments_{timestamp}.csv"
#             # df.to_csv(filename, index=False, encoding='utf-8')
#             # print(f"✅ Данные сохранены в: {filename}")
#             # print(f"📊 Статистика:")
#             # print(f"   - Всего фрагментов: {len(df)}")
#             # print(f"   - Диапазон размеров: {start_size}-{end_size} нт")
#         # return pd.DataFrame(self._gc_content)
#
#
#
#     # # Создаем DataFrame с результатами
#     #
#     #
#     # # Выводим красивые результаты
#     # print("Результаты анализа GC-content:")
#     # print(df[['sequence_id', 'length', 'gc_percent', 'gc_category']].to_string(index=False))
#     #
#     # # Статистика
#     # print(f"\nСтатистика:")
#     # print(f"Средний GC-content: {df['gc_percent'].mean():.1f}%")
#     # print(f"Минимальный GC-content: {df['gc_percent'].min():.1f}%")
#     # print(f"Максимальный GC-content: {df['gc_percent'].max():.1f}%")
#     #
#     # with pd.ExcelWriter('gc_analys.xlsx') as writer:
#     #     # Основные данные
#     #     df.to_excel(writer, sheet_name='Sequences', index=False)
#     #
#     #     # Статистика
#     #     df.to_excel(writer, sheet_name='Statistics', index=False)
#     #
#     # print("Файл 'gc_analys.xlsx' успешно создан!")


