# def range_sequence(self):
    #     """Функция обрезает последовательность РНК по задаваемому диапазону"""
    #     # start = int(input("Введите начальный индекс: "))
    #     # end = int(input("Введите конечный индекс: "))
    #     # seq1 = self.sequence[start:end]
    #     # print(f"Срез от {start} до {end}:")
    #     # print(seq1)
    #     # print(f"Длина среза: {len(seq1)}")
    #      = self.sequence[781:2858]
    #     return


        # df = pd.DataFrame(fragments_list)
        #
        # df.to_csv(filename, index=False, encoding='utf-8')
        #
        # print(f"✅ Данные сохранены в: {filename}")
        # print(f"📊 Статистика:")
        # print(f"   - Всего фрагментов: {len(df)}")
        # print(f"   - Диапазон размеров: {start_size}-{end_size} нт")
        # print(f"   - Размеры в файле: {df['size_nt':fra].unique()}")

def smart_merge_csv(files_list, output_file):
    """
    Умное объединение по общим столбцам
    """
    dataframes = []

    for file_path in files_list:
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            dataframes.append(df)
            print(f"Загружен {file_path}: {df.shape}")
        else:
            print(f"Файл {file_path} не найден")

    if not dataframes:
        print("Нет файлов для объединения")
        return

    # Находим общие столбцы
    common_columns = set(dataframes[0].columns)
    for df in dataframes[1:]:
        common_columns = common_columns.intersection(set(df.columns))

    print(f"Общие столбцы: {list(common_columns)}")

    if not common_columns:
        print("Нет общих столбцов для объединения")
        return

    # Объединяем по общим столбцам
    merged_df = dataframes[0]
    for df in dataframes[1:]:
        merged_df = pd.merge(merged_df, df, on=list(common_columns), how='outer')

    merged_df.to_csv(output_file, index=False)
    print(f"Объединенный файл сохранен: {output_file}")
    return merged_df