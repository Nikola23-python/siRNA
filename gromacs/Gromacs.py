import os
import subprocess
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import shutil

PDB_DIR = "all_pdbs"
OUTPUT_FILE = "final_energies.csv"
GMX = "gmx"

TEMP_DIR = "temp_gromacs"
Path(TEMP_DIR).mkdir(exist_ok=True)


def run_minimization(pdb_file):
    work_dir = Path(TEMP_DIR) / pdb_file.stem
    work_dir.mkdir(exist_ok=True)
    shutil.copy(pdb_file, work_dir / "input.pdb")
    original_dir = os.getcwd()
    os.chdir(work_dir)

    try:
        cmd_topo = f"echo 1 1 | {GMX} pdb2gmx -f input.pdb -o conf.gro -ff amber99sb-ildn -water tip3p -ignh -ter"
        subprocess.run(cmd_topo, shell=True, capture_output=True, text=True)
        subprocess.run(f"{GMX} editconf -f conf.gro -o box.gro -c -d 1.0",
                       shell=True, capture_output=True)
        subprocess.run(f"{GMX} solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top",
                       shell=True, capture_output=True)
        with open("minim.mdp", "w") as f:
            f.write("integrator=steep\n")
            f.write("emtol=1000.0\n")
            f.write("nsteps=5000\n")
            f.write("cutoff-scheme=Verlet\n")
            f.write("coulombtype=PME\n")
            f.write("rcoulomb=1.0\n")
            f.write("rvdw=1.0\n")
            f.write("pbc=xyz\n")

        subprocess.run(f"{GMX} grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 5",
                       shell=True, capture_output=True)
        subprocess.run(f"{GMX} mdrun -deffnm em -nt 1",
                       shell=True, capture_output=True)

        subprocess.run(f'echo "Potential" | {GMX} energy -f em.edr -o en.xvg',
                       shell=True, capture_output=True)

        energy = None
        if Path("en.xvg").exists():
            with open("en.xvg") as f:
                for line in f:
                    if not line.startswith(('@', '#')):
                        energy = float(line.split()[1])

        return energy

    except Exception as e:
        print(f"   Ошибка: {e}")
        return None
    finally:
        os.chdir(original_dir)


def parse_filename(filename):
    name = filename.stem
    parts = name.split('_')
    for i, part in enumerate(parts):
        if part == "si":
            idx = parts[i + 1]
            oh_type = parts[i + 2]
            return idx, oh_type
    return "unknown", "unknown"


def main():
    pdb_files = list(Path(PDB_DIR).glob("*.pdb"))
    if not pdb_files:
        print(f"❌ Нет PDB файлов в папке {PDB_DIR}/")
        return
    print(f" Найдено PDB файлов: {len(pdb_files)}")

    results = []

    for pdb_file in tqdm(pdb_files, desc="Минимизация"):
        idx, oh_type = parse_filename(pdb_file)

        energy = run_minimization(pdb_file)

        results.append({
            "id": idx,
            "overhang": oh_type,
            "energy": energy,
            "file": pdb_file.name
        })
        time.sleep(1)

    df = pd.DataFrame(results)
    df.to_csv(OUTPUT_FILE, index=False)

    success = df['energy'].notna().sum()
    print(f"\n📊 Результаты:")
    print(f"   Всего: {len(df)}")
    print(f"   Успешно: {success}")
    print(f"   Ошибок: {len(df) - success}")
    print(f"\n💾 Результаты сохранены в {OUTPUT_FILE}")

if __name__ == "__main__":
    import time
    main()