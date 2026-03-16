import os
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import subprocess

# --- НАСТРОЙКИ ---
PDB_DIR = "all_pdbs"
EXISTING_FILE = "final_energies_clean.csv"
OUTPUT_FILE = "final_energies_complete.csv"
GMX = "gmx"
TEMP_DIR = "temp_gromacs"

Path(TEMP_DIR).mkdir(exist_ok=True)

# Загружаем уже готовые
df_existing = pd.read_csv(EXISTING_FILE)
existing_set = set(zip(df_existing['id'], df_existing['overhang']))
print(f"✅ Уже есть: {len(existing_set)}")

# Все PDB файлы
pdb_files = list(Path(PDB_DIR).glob("*.pdb"))
print(f"📁 Всего PDB: {len(pdb_files)}")


def parse_filename(filename):
    name = filename.stem
    parts = name.split('_')
    for i, part in enumerate(parts):
        if part == "si":
            try:
                idx = int(parts[i + 1])
                oh_type = parts[i + 2]
                return idx, oh_type
            except:
                pass
    return None, None


def run_minimization(pdb_path, id_num, oh_type):
    work_dir = Path(TEMP_DIR) / f"si_{id_num}_{oh_type}"
    work_dir.mkdir(exist_ok=True)

    os.system(f"cp '{pdb_path}' '{work_dir}/input.pdb' 2>/dev/null")
    original_dir = os.getcwd()
    os.chdir(work_dir)

    try:
        # pdb2gmx
        subprocess.run(
            f"echo '1 1' | {GMX} pdb2gmx -f input.pdb -o conf.gro -ff amber99sb-ildn -water tip3p -ignh -ter",
            shell=True, capture_output=True)

        # editconf
        subprocess.run(f"{GMX} editconf -f conf.gro -o box.gro -c -d 1.0",
                       shell=True, capture_output=True)

        # solvate
        subprocess.run(f"{GMX} solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top",
                       shell=True, capture_output=True)

        # minim.mdp
        with open("minim.mdp", "w") as f:
            f.write("integrator=steep\n")
            f.write("emtol=1000.0\n")
            f.write("nsteps=5000\n")
            f.write("cutoff-scheme=Verlet\n")

        # grompp
        subprocess.run(f"{GMX} grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 5",
                       shell=True, capture_output=True)

        # mdrun
        subprocess.run(f"{GMX} mdrun -deffnm em -nt 1",
                       shell=True, capture_output=True)

        # energy
        subprocess.run(f'echo "Potential" | {GMX} energy -f em.edr -o en.xvg',
                       shell=True, capture_output=True)

        # парсим
        if Path("en.xvg").exists():
            with open("en.xvg") as f:
                for line in f:
                    if not line.startswith(('@', '#')) and line.strip():
                        energy = float(line.split()[1])
                        os.chdir(original_dir)
                        return energy
    except Exception as e:
        print(f"   Ошибка: {e}")

    os.chdir(original_dir)
    return None


# Собираем только недостающие
to_process = []
for pdb_file in pdb_files:
    id_num, oh_type = parse_filename(pdb_file)
    if id_num is not None and (id_num, oh_type) not in existing_set:
        to_process.append((pdb_file, id_num, oh_type))

print(f"🎯 Нужно обработать: {len(to_process)}")

results = []
for pdb_file, id_num, oh_type in tqdm(to_process, desc="Минимизация"):
    energy = run_minimization(pdb_file, id_num, oh_type)
    if energy is not None:
        results.append({"id": id_num, "overhang": oh_type, "energy": energy})

    # Сохраняем каждые 20 файлов
    if len(results) % 20 == 0 and results:
        df_temp = pd.DataFrame(results)
        df_all = pd.concat([df_existing, df_temp], ignore_index=True)
        df_all.to_csv(OUTPUT_FILE, index=False)
        print(f"\n💾 Промежуточно сохранено: {len(df_all)}")

# Финальное сохранение
df_new = pd.DataFrame(results)
df_final = pd.concat([df_existing, df_new], ignore_index=True)
df_final = df_final.sort_values(['id', 'overhang']).reset_index(drop=True)
df_final.to_csv(OUTPUT_FILE, index=False)

print(f"\n{'=' * 50}")
print(f"✅ Всего в датасете: {len(df_final)}")
print(f"📊 Распределение по оверхенгам:")
print(df_final['overhang'].value_counts())
print(f"{'=' * 50}")
print(f"💾 Финальный файл: {OUTPUT_FILE}")