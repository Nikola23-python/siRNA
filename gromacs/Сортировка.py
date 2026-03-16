import zipfile
import os
import shutil
from pathlib import Path

# --- НАСТРОЙКИ ---
RESULTS_DIR = "rna_compozer"  # Папка куда скачиваются архивы
OUTPUT_DIR = "all_pdbs"  # Папка куда сложим все PDB
PROCESSED_DIR = "processed"  # Папка для обработанных архивов

# Создаем папки
Path(OUTPUT_DIR).mkdir(exist_ok=True)
Path(PROCESSED_DIR).mkdir(exist_ok=True)


def extract_pdbs_from_zip(zip_path):
    temp_dir = Path("temp_extract")
    temp_dir.mkdir(exist_ok=True)

    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(temp_dir)
            print(f"   Распаковано {len(zip_ref.namelist())} файлов")
        pdb_files = list(temp_dir.rglob("*.pdb"))
        print(f"   Найдено PDB: {len(pdb_files)}")

        for pdb_file in pdb_files:
            new_name = f"{zip_path.stem}_{pdb_file.name}"
            dest_path = Path(OUTPUT_DIR) / new_name
            shutil.copy2(pdb_file, dest_path)
            print(f"   ➜ {new_name}")
        shutil.move(str(zip_path), Path(PROCESSED_DIR) / zip_path.name)
        print(f"   ✅ Архив перемещен в {PROCESSED_DIR}/")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


def main():
    zip_files = list(Path(RESULTS_DIR).glob("*.zip"))
    if not zip_files:
        print("❌ Нет zip файлов в папке results/")
        return
    print(f"Найдено архивов: {len(zip_files)}")
    for zip_file in zip_files:
        extract_pdbs_from_zip(zip_file)
    all_pdbs = list(Path(OUTPUT_DIR).glob("*.pdb"))
    print(f"\nВсего собрано PDB файлов: {len(all_pdbs)}")
    print(f"Все PDB в папке: {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()