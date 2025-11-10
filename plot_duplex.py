import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import pandas as pd
from prepare_rna import df_sense, df_antisense
from rna_duplex import duplex_df


def get_sequences_by_id(sense_id, antisense_id):
    """Получить последовательности по ID"""
    sense_seq = df_sense[df_sense['fragment_id'] == sense_id]['sequence'].iloc[0]
    antisense_seq = df_antisense[df_antisense['fragment_id'] == antisense_id]['sequence'].iloc[0]
    return sense_seq, antisense_seq


def get_duplex_data(sense_id, antisense_id):
    """Получить данные дуплекса по ID"""
    duplex_data = duplex_df[
        (duplex_df['sense_id'] == sense_id) &
        (duplex_df['antisense_id'] == antisense_id)
        ].iloc[0]
    return duplex_data


def get_base_color(base, symbol):
    """Возвращает цвет для основания"""
    base_colors = {
        'A': '#FF6B6B', 'U': '#4ECDC4', 'G': '#FFD166', 'C': '#06D6A0'
    }

    if symbol == '(':
        return '#FF5252'
    elif symbol == ')':
        return '#448AFF'
    else:
        return base_colors.get(base, '#CFCFCF')


def calculate_connections(sense_seq, antisense_seq, structure):
    """Вычисляет связи между цепями"""
    connections = []
    stack = []
    sense_length = len(sense_seq)

    for i, symbol in enumerate(structure):
        if symbol == '&':
            continue
        elif symbol == '(':
            stack.append(i)
        elif symbol == ')':
            if stack:
                j = stack.pop()
                if j < sense_length and i >= sense_length + 1:
                    sense_pos = j
                    anti_pos = i - sense_length - 1
                    if anti_pos < len(antisense_seq):
                        connections.append((sense_pos, anti_pos))
    return connections


def plot_advanced_structure(ax, sense_seq, antisense_seq, structure, energy, title):
    """Продвинутая визуализация структуры с красивым оформлением"""
    parts = structure.split('&')
    sense_struct = parts[0]
    antisense_struct = parts[1]

    ax.clear()
    max_len = max(len(sense_seq), len(antisense_seq))
    base_spacing = 1.0
    vertical_spacing = 2.0

    y_sense = vertical_spacing
    sense_start_x = (max_len * base_spacing) / 2 - (len(sense_seq) * base_spacing) / 2
    y_antisense = 0
    antisense_start_x = (max_len * base_spacing) / 2 - (len(antisense_seq) * base_spacing) / 2

    connections = calculate_connections(sense_seq, antisense_seq, structure)
    for sense_pos, anti_pos in connections:
        x1 = sense_start_x + sense_pos * base_spacing
        x2 = antisense_start_x + anti_pos * base_spacing
        ax.plot([x1, x2], [y_sense - 0.15, y_antisense + 0.15],
                color='purple', alpha=0.7, linewidth=2.5, zorder=1)

    for i, (base, symbol) in enumerate(zip(sense_seq, sense_struct)):
        x_pos = sense_start_x + i * base_spacing
        color = get_base_color(base, symbol)
        patch = FancyBboxPatch((x_pos - 0.3, y_sense - 0.3), 0.6, 0.6,
                               boxstyle="round,pad=0.02",
                               facecolor=color, edgecolor='black',
                               linewidth=1.5, zorder=2)
        ax.add_patch(patch)
        ax.text(x_pos, y_sense, base, ha='center', va='center',
                fontsize=14, fontweight='bold', zorder=3)

    for i, (base, symbol) in enumerate(zip(antisense_seq, antisense_struct)):
        x_pos = antisense_start_x + i * base_spacing
        color = get_base_color(base, symbol)
        patch = FancyBboxPatch((x_pos - 0.3, y_antisense - 0.3), 0.6, 0.6,
                               boxstyle="round,pad=0.02",
                               facecolor=color, edgecolor='black',
                               linewidth=1.5, zorder=2)
        ax.add_patch(patch)
        ax.text(x_pos, y_antisense, base, ha='center', va='center',
                fontsize=14, fontweight='bold', zorder=3)

    ax.text(sense_start_x - 1, y_sense, 'SENSE', ha='right', va='center',
            fontsize=12, fontweight='bold', color='darkred')
    ax.text(antisense_start_x - 1, y_antisense, 'ANTI', ha='right', va='center',
            fontsize=12, fontweight='bold', color='darkblue')
    ax.set_title(f'{title}\nEnergy: {energy:.2f} kcal/mol',
                 fontsize=16, fontweight='bold', pad=20)
    margin = 2
    ax.set_xlim(-margin, max_len * base_spacing + margin)
    ax.set_ylim(-1, vertical_spacing + 1)
    ax.set_aspect('equal')
    ax.axis('off')


def plot_enhanced_info(ax, sense_seq, antisense_seq, structure, energy, sense_id, antisense_id):
    """Улучшенная информационная панель"""
    ax.clear()
    ax.axis('off')
    gc_sense = (sense_seq.count('G') + sense_seq.count('C')) / len(sense_seq) * 100
    gc_anti = (antisense_seq.count('G') + antisense_seq.count('C')) / len(antisense_seq) * 100
    stability = "POOR" if energy > -5 else "MODERATE" if energy > -10 else "GOOD" if energy > -15 else "EXCELLENT"

    info_text = f"""
    DUPLEX ANALYSIS REPORT
    ----------------------
    IDs:
    Sense: {sense_id}
    Antisense: {antisense_id}

    Sequences:
    Sense:    {sense_seq}
    Antisense: {antisense_seq}

    STRUCTURE: 
    {structure}

    STATISTICS:
    Sense length: {len(sense_seq)} nt
    Anti length:  {len(antisense_seq)} nt
    GC Sense: {gc_sense:.1f}%
    GC Anti:  {gc_anti:.1f}%
    Stability: {stability}
    Energy: {energy:.2f} kcal/mol
    Base pairs: {len(calculate_connections(sense_seq, antisense_seq, structure))}
    """

    ax.text(0.02, 0.98, info_text, transform=ax.transAxes, fontsize=12,
            fontfamily='DejaVu Sans Mono', verticalalignment='top', linespacing=1.5,
            bbox=dict(boxstyle="round,pad=1", facecolor='lightblue', alpha=0.2))


def plot_duplex_by_id(sense_id, antisense_id, save_path=None, figsize=(18, 8)):
    """График дуплекса по ID"""
    sense_seq, antisense_seq = get_sequences_by_id(sense_id, antisense_id)
    duplex_data = get_duplex_data(sense_id, antisense_id)
    structure = duplex_data['duplex_structure']
    energy = duplex_data['duplex_energy']

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
    ax_struct = plt.subplot(gs[0])
    ax_info = plt.subplot(gs[1])

    plot_advanced_structure(ax_struct, sense_seq, antisense_seq, structure, energy, f"{sense_id} + {antisense_id}")
    plot_enhanced_info(ax_info, sense_seq, antisense_seq, structure, energy, sense_id, antisense_id)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"График сохранен: {save_path}")
    plt.show()
    return duplex_data


def show_available_sequences():
    """Показать доступные последовательности"""
    print("\nДоступные sense последовательности:")
    sense_ids = df_sense['fragment_id'].unique()
    for i, sid in enumerate(sense_ids[:20]):  # Показываем первые 20
        print(f"  {sid}")
    if len(sense_ids) > 20:
        print(f"  ... и еще {len(sense_ids) - 20} последовательностей")

    print("\nДоступные antisense последовательности:")
    anti_ids = df_antisense['fragment_id'].unique()
    for i, aid in enumerate(anti_ids[:20]):  # Показываем первые 20
        print(f"  {aid}")
    if len(anti_ids) > 20:
        print(f"  ... и еще {len(anti_ids) - 20} последовательностей")


def interactive_plotter():
    """Интерактивный построитель графиков"""
    print("=== ИНТЕРАКТИВНЫЙ ПОСТРОИТЕЛЬ ГРАФИКОВ DUPLEX ===")

    while True:
        print("\n" + "=" * 50)
        print("1. Показать доступные последовательности")
        print("2. Построить график дуплекса")
        print("3. Выход")

        choice = input("\nВыберите опцию (1-3): ").strip()

        if choice == '1':
            show_available_sequences()

        elif choice == '2':
            print("\nВведите ID последовательностей:")
            sense_id = input("Sense ID: ").strip()
            antisense_id = input("Antisense ID: ").strip()

            # Проверяем существование последовательностей
            try:
                sense_seq, antisense_seq = get_sequences_by_id(sense_id, antisense_id)
                print(f"Найдены последовательности:")
                print(f"  Sense: {sense_seq}")
                print(f"  Antisense: {antisense_seq}")

                # Запрос имени файла для сохранения
                filename = input("\nВведите имя файла для сохранения (или нажмите Enter для пропуска): ").strip()
                save_path = None
                if filename:
                    if not filename.endswith('.png'):
                        filename += '.png'
                    save_path = filename
                    print(f"График будет сохранен как: {save_path}")

                # Построение графика
                print("\nСтроим график...")
                plot_duplex_by_id(sense_id, antisense_id, save_path)

            except Exception as e:
                print(f"Ошибка: {e}")
                print("Проверьте правильность ID последовательностей")

        elif choice == '3':
            print("Выход из программы")
            break

        else:
            print("Неверный выбор. Попробуйте снова.")


if __name__ == "__main__":
    interactive_plotter()