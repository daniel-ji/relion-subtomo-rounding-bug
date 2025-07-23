"""
Standalone script to compare two MRC files for subtomogram extraction consistency.
"""

import argparse
import mrcfile
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def print_statistics(name: str, data: np.ndarray) -> None:
    median = np.median(data)
    p25, p75 = np.percentile(data, [25, 75])
    print(f"  - {name}: Median={median:.4f}, IQR=({p25:.4f} to {p75:.4f})")

def plot_heatmaps(
    mrc_file_1_name: str,
    mrc_file_2_name: str,
    mrc_file_1_data: np.ndarray,
    mrc_file_2_data: np.ndarray,
    difference: np.ndarray,
    section_num: int,
    output_dir: Path,
) -> None:
    heatmap_kwargs = {"xticklabels": False, "yticklabels": False}
    vmin_1, vmax_1 = np.percentile(mrc_file_1_data, [5, 95])
    vmin_2, vmax_2 = np.percentile(mrc_file_2_data, [5, 95])
    vmin_diff, vmax_diff = np.percentile(difference, [5, 95])

    fig, axes = plt.subplots(figsize=(10, 6))
    sns.heatmap(mrc_file_1_data, ax=axes, cmap='viridis', vmin=vmin_1, vmax=vmax_1, **heatmap_kwargs).set_title(f'{mrc_file_1_name} Data')
    plt.savefig(output_dir / f'section_{section_num}_{mrc_file_1_name}_data.png')
    plt.close(fig)

    fig, axes = plt.subplots(figsize=(10, 6))
    sns.heatmap(mrc_file_2_data, ax=axes, cmap='viridis', vmin=vmin_2, vmax=vmax_2, **heatmap_kwargs).set_title(f'{mrc_file_2_name} Data')
    plt.savefig(output_dir / f'section_{section_num}_{mrc_file_2_name}_data.png')
    plt.close(fig)

    fig, axes = plt.subplots(figsize=(10, 6))
    sns.heatmap(difference, ax=axes, cmap='hot', vmin=vmin_diff, vmax=vmax_diff, **heatmap_kwargs).set_title('Absolute Difference')
    plt.savefig(output_dir / f'section_{section_num}_absolute_difference.png')
    plt.close(fig)


def compare(
    data1_filename: str,
    data2_filename: str,
    data1: np.ndarray,
    data2: np.ndarray,
    section_num: int,
    output_dir: Path,
) -> None:
    difference = np.abs(data1 - data2)

    print_statistics(f"{data1_filename} Data", data1)
    print_statistics(f"{data2_filename} Data", data2)
    print_statistics("Difference", difference)
    plot_heatmaps(data1_filename, data2_filename, data1, data2, difference, section_num, output_dir)

def main():
    parser = argparse.ArgumentParser(
        description="Compare 2D sections of two MRC files.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--mrc-file-1-name",
        type=str,
        required=True,
        help="Name of the first MRC file for display purposes.",
    )
    parser.add_argument(
        "--mrc-file-2-name",
        type=str,
        required=True,
        help="Name of the second MRC file for display purposes.",
    )
    parser.add_argument(
        "--mrc-file-1",
        type=Path,
        required=True,
        help="Path to the first MRC file.",
    )
    parser.add_argument(
        "--mrc-file-2",
        type=Path,
        required=True,
        help="Path to the second MRC file.",
    )
    parser.add_argument(
        "--sections",
        type=int,
        nargs='+',
        required=True,
        help="One or more space-separated section numbers to compare (e.g., 0 5 10).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("./mrc_comparison_results"),
        help="Directory to save the output plots.",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    with mrcfile.open(args.mrc_file_1, permissive=True) as mrc_file_1, \
            mrcfile.open(args.mrc_file_2, permissive=True) as mrc_file_2:

        mrc_file_1_data = mrc_file_1.data
        mrc_file_2_data = mrc_file_2.data

        if mrc_file_1_data.shape != mrc_file_2_data.shape:
            raise ValueError(
                f"MRC files must have the same shape. "
                f"File 1: {mrc_file_1_data.shape}, File 2: {mrc_file_2_data.shape}"
            )

        for section in args.sections:
            print(f"\n{'='*20} Processing Section {section} {'='*20}")
            if section > mrc_file_1_data.shape[0]:
                print(f"Warning: Section {section} is out of bounds for shape {mrc_file_1_data.shape}. Skipping.")
                continue
            
            compare(
                args.mrc_file_1_name,
                args.mrc_file_2_name,
                mrc_file_1_data[section - 1],
                mrc_file_2_data[section - 1],
                section,
                args.output_dir,
            )

if __name__ == "__main__":
    main()