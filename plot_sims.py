#!/usr/bin/env python3
"""Render vis_v2's linear and log-count PDFs from simulation CSVs.

The plotting code is copied from vis_v2.ipynb, including its styles, axis
limits, interpolation, filenames, and newest-file-per-scenario selection.
Use --animations to also export the notebook's GIFs. No notebook is needed.
Requires numpy, pandas, matplotlib, and a working LaTeX installation;
animations additionally use Pillow (normally installed with matplotlib).
"""

import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir", type=Path, default=Path("final_submission_sims/hill_sims"),
        help="Directory containing simulation CSVs (default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        help="Directory for PDFs and optional GIFs (default: DATA_DIR/figures).",
    )
    parser.add_argument(
        "--animations", action="store_true",
        help="Also create the notebook's linear and log-count animated GIFs.",
    )
    args = parser.parse_args()
    if not args.data_dir.is_dir():
        parser.error(f"Input directory does not exist: {args.data_dir}")
    if not any(args.data_dir.glob("simulation_history_*.csv")):
        parser.error(f"No simulation CSVs found in {args.data_dir}")

    import matplotlib
    matplotlib.use("Agg")

    import re

    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    DATA_DIR = args.data_dir
    FIGURES_DIR = args.output_dir if args.output_dir is not None else DATA_DIR / "figures"
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Arial"],
        "mathtext.fontset": "cm",
        "font.size": 12,
        "axes.labelsize": 12,
        "legend.fontsize": 11,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.dpi": 200,
        "axes.titlesize": 9,
        "axes.linewidth": 0.6,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "legend.frameon": False,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 2,
        "ytick.minor.size": 2,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "grid.alpha": 0.3,
        "savefig.dpi": 300,
        "savefig.format": "pdf",
        "savefig.bbox": "tight",
    })

    PATTERNS = [
        (
            "regular",
            re.compile(
                r"^simulation_history_n(\d+)"
                r"_beta([\d.]+)"
                r"_gamma([\d.]+)"
                r"_rhoT([\d.]+)"
                r"_rhoH([\d.]+)"
                r"_c([\d.]+)"
                r"_maxTime([\d.]+)"
                r"(?:_seed(\d+))?"
                r"\.csv$"
            ),
            "",
        ),
        (
            "travel_ban",
            re.compile(
                r"^simulation_history_n(\d+)"
                r"_beta([\d.]+)"
                r"_gamma([\d.]+)"
                r"_rhoT([\d.]+)"
                r"_rhoH([\d.]+)"
                r"_c([\d.]+)"
                r"_maxTime([\d.]+)"
                r"(?:_seed(\d+))?"
                r"_TB"
                r"\.csv$"
            ),
            "_TB",
        ),
        (
            "social_distancing",
            re.compile(
                r"^simulation_history_n(\d+)"
                r"_beta([\d.]+)"
                r"_gamma([\d.]+)"
                r"_rhoT([\d.]+)"
                r"_rhoH([\d.]+)"
                r"_c([\d.]+)"
                r"_maxTime([\d.]+)"
                r"(?:_seed(\d+))?"
                r"_SD([\d.]+)"
                r"\.csv$"
            ),
            "_SD",
        ),
    ]

    x_upper_lim = None

    for label, pattern, suffix in PATTERNS:
        candidates = [path for path in DATA_DIR.iterdir() if pattern.match(path.name)]
        if not candidates:
            print(f"Skipping {label}: no matching files in {DATA_DIR}")
            continue

        latest_file = max(candidates, key=lambda path: path.stat().st_mtime)
        match = pattern.match(latest_file.name)

        n = int(match.group(1))
        beta = float(match.group(2))
        gamma = float(match.group(3))
        rho_T = float(match.group(4))
        rho_H = float(match.group(5))
        c = int(match.group(6))
        max_time = int(match.group(7))
        seed = match.group(8) or "0"

        print(f"Reading {label}: {latest_file}")

        df = pd.read_csv(latest_file)
        end_rows = df[(df["I1"] == 0) & (df["I2"] == 0)].index
        if len(end_rows) > 0:
            df = df.iloc[: end_rows[0] + 5]

        time_vals = df["time"].values
        I1 = df["I1"].values
        R1 = df["R1"].values
        I2 = df["I2"].values
        R2 = df["R2"].values

        tau_scale = np.log(n)
        adjusted_times = time_vals / tau_scale

        if x_upper_lim is None:
            x_upper_lim = max(adjusted_times)
        elif x_upper_lim > max(adjusted_times):
            adjusted_times = np.append(adjusted_times, x_upper_lim)
            I1 = np.append(I1, I1[-1])
            I2 = np.append(I2, I2[-1])
            R1 = np.append(R1, R1[-1])
            R2 = np.append(R2, R2[-1])

        fig, ax1 = plt.subplots(figsize=(10, 6))

        ax1.plot(adjusted_times, I1 / n * 100, label="Infected (Community 1)", color="C0", linestyle="--")
        ax1.plot(adjusted_times, I2 / n * 100, label="Infected (Community 2)", color="C1", linestyle="--")
        ax1.set_xlabel(r"Scaled time ($\ln n$)")
        ax1.set_ylabel(r"Percentage of population (\%)")
        ax1.set_xlim(0, x_upper_lim)
        ax1.set_ylim(0, 30)

        ax2 = ax1.twinx()
        ax2.fill_between(adjusted_times, 0, (I1 + R1) / n * 100, color="C0", alpha=0.1)
        ax2.fill_between(adjusted_times, 0, (I2 + R2) / n * 100, color="C1", alpha=0.1)
        ax2.plot(adjusted_times, (I1 + R1) / n * 100, color="C0", linestyle="-", linewidth=0.1)
        ax2.plot(adjusted_times, (I2 + R2) / n * 100, color="C1", linestyle="-", linewidth=0.1)
        ax2.set_ylabel(r"Cumulative percentage infected (\%)")
        ax2.set_ylim(0, 30)

        ax1.legend(loc="upper right")
        fig.subplots_adjust(left=0.12, right=0.88, bottom=0.12, top=0.96)

        output_file = FIGURES_DIR / (
            f"simulation_n{n}_b{beta}_g{gamma}_T{rho_T}_H{rho_H}_c{c}_t{max_time}_seed{seed}{suffix}.pdf"
        )
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved {output_file}")
        plt.close(fig)

    def load_latest_simulation(pattern):
        candidates = [path for path in DATA_DIR.iterdir() if pattern.match(path.name)]
        if not candidates:
            return None

        latest_file = max(candidates, key=lambda path: path.stat().st_mtime)
        match = pattern.match(latest_file.name)

        n = int(match.group(1))
        beta = float(match.group(2))
        gamma = float(match.group(3))
        rho_T = float(match.group(4))
        rho_H = float(match.group(5))
        c = int(match.group(6))
        max_time = int(match.group(7))
        seed = match.group(8) or "0"

        df = pd.read_csv(latest_file)
        end_rows = df[(df["I1"] == 0) & (df["I2"] == 0)].index
        if len(end_rows) > 0:
            df = df.iloc[: end_rows[0] + 5]

        return {
            "file": latest_file,
            "n": n,
            "beta": beta,
            "gamma": gamma,
            "rho_T": rho_T,
            "rho_H": rho_H,
            "c": c,
            "max_time": max_time,
            "seed": seed,
            "time": df["time"].values / np.log(n),
            "I1": df["I1"].values,
            "I2": df["I2"].values,
            "R1": df["R1"].values,
            "R2": df["R2"].values,
        }

    records = []
    for label, pattern, suffix in PATTERNS:
        record = load_latest_simulation(pattern)
        if record is None:
            print(f"Skipping {label}: no matching files in {DATA_DIR}")
            continue

        record["label"] = label
        record["suffix"] = suffix
        records.append(record)

    if not records:
        raise ValueError(f"No simulation files found in {DATA_DIR}")

    x_upper_lim = max(record["time"].max() for record in records)
    y_upper_lim = 1e8

    for record in records:
        print(f"Reading {record['label']}: {record['file']}")

        time = record["time"].copy()
        I1 = record["I1"].copy()
        I2 = record["I2"].copy()
        C1 = (record["I1"] + record["R1"]).copy()
        C2 = (record["I2"] + record["R2"]).copy()

        # Extend shorter simulations to the shared x-axis limit so the cumulative fill is not cut off.
        if time[-1] < x_upper_lim:
            time = np.append(time, x_upper_lim)
            I1 = np.append(I1, I1[-1])
            I2 = np.append(I2, I2[-1])
            C1 = np.append(C1, C1[-1])
            C2 = np.append(C2, C2[-1])

        I1_log = np.where(I1 > 0, I1, np.nan)
        I2_log = np.where(I2 > 0, I2, np.nan)
        C1_log = np.where(C1 > 0, C1, np.nan)
        C2_log = np.where(C2 > 0, C2, np.nan)

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.fill_between(time, 1, np.maximum(C1, 1), where=C1 > 0, color="C0", alpha=0.10, linewidth=0)
        ax.fill_between(time, 1, np.maximum(C2, 1), where=C2 > 0, color="C1", alpha=0.10, linewidth=0)

        ax.plot(time, I1_log, label="Infected (Community 1)", color="C0", linestyle="--", linewidth=1)
        ax.plot(time, I2_log, label="Infected (Community 2)", color="C1", linestyle="--", linewidth=1)

        ax.plot(time, C1_log, color="C0", linestyle="-", linewidth=0.1)
        ax.plot(time, C2_log, color="C1", linestyle="-", linewidth=0.1)

        # Mark the first seed into Community 2 on the intervention plots only.
        first_seed_indices = np.flatnonzero(C2 > 0)
        if record["suffix"] and len(first_seed_indices) > 0:
            first_seed_time = time[first_seed_indices[0]]
            ax.axvline(
                first_seed_time,
                color="0.45",
                linestyle="--",
                linewidth=1.2,
                label="First seed",
            )

        threshold = 0.01 * record["n"]
        threshold_indices = np.flatnonzero(C1 >= threshold)
        if record["suffix"] and len(threshold_indices) > 0:
            idx = threshold_indices[0]
            if idx == 0 or C1[idx] == C1[idx - 1]:
                threshold_time = time[idx]
            else:
                threshold_time = time[idx - 1] + (
                    (threshold - C1[idx - 1])
                    * (time[idx] - time[idx - 1])
                    / (C1[idx] - C1[idx - 1])
                )
            ax.axvline(threshold_time, color="red", linestyle="--", linewidth=1.0, label=r"$\tau_1(1\%)$")

        ax.set_yscale("log")
        ax.set_xlabel(r"Scaled time ($\ln n$)")
        ax.set_ylabel("Number of individuals")
        ax.set_xlim(0, x_upper_lim)
        ax.set_ylim(1, y_upper_lim * 1.2)
        ax.legend(loc="upper right")

        fig.subplots_adjust(left=0.14, right=0.96, bottom=0.12, top=0.96)

        output_file = FIGURES_DIR / (
            f"simulation_n{record['n']}_b{record['beta']}_g{record['gamma']}_T{record['rho_T']}_H{record['rho_H']}_"
            f"c{record['c']}_t{record['max_time']}_seed{record['seed']}{record['suffix']}_log_counts.pdf"
        )
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        print(f"Saved {output_file}")
        plt.close(fig)

    if args.animations:
        from matplotlib.animation import FuncAnimation, PillowWriter

        ANIMATION_FPS = 16
        ANIMATION_FRAMES = 120

        SCENARIO_TITLES = {
            "regular": "No intervention",
            "travel_ban": "Travel ban",
            "social_distancing": "Social distancing",
        }

        def load_animation_records():
            records = []
            for label, pattern, suffix in PATTERNS:
                record = load_latest_simulation(pattern)
                if record is None:
                    print(f"Skipping {label}: no matching files in {DATA_DIR}")
                    continue
                df = pd.read_csv(record["file"])
                end_rows = df[(df["I1"] == 0) & (df["I2"] == 0)].index
                if len(end_rows) > 0:
                    df = df.iloc[: end_rows[0] + 5]
                record.update({
                    "label": label,
                    "suffix": suffix,
                    "title": SCENARIO_TITLES.get(label, label.replace("_", " ").title()),
                    "time": df["time"].values / np.log(record["n"]),
                    "I1": df["I1"].values,
                    "I2": df["I2"].values,
                    "R1": df["R1"].values,
                    "R2": df["R2"].values,
                })
                records.append(record)
            if not records:
                raise ValueError(f"No simulation files found in {DATA_DIR}")
            return records

        def first_threshold_time(time, infected, n):
            threshold = 0.01 * n
            threshold_indices = np.flatnonzero(infected >= threshold)
            if len(threshold_indices) == 0:
                return None
            idx = threshold_indices[0]
            if idx == 0 or infected[idx] == infected[idx - 1]:
                return time[idx]
            return time[idx - 1] + (threshold - infected[idx - 1]) * (time[idx] - time[idx - 1]) / (infected[idx] - infected[idx - 1])

        animation_records = load_animation_records()
        animation_x_upper_lim = max(record["time"].max() for record in animation_records)

        for record in animation_records:
            time = record["time"]
            I1_pct = record["I1"] / record["n"] * 100
            I2_pct = record["I2"] / record["n"] * 100
            C1_pct = (record["I1"] + record["R1"]) / record["n"] * 100
            C2_pct = (record["I2"] + record["R2"]) / record["n"] * 100
            frame_indices = np.unique(np.linspace(2, len(time), min(ANIMATION_FRAMES, len(time)), dtype=int))

            fig, ax1 = plt.subplots(figsize=(10, 6))
            ax2 = ax1.twinx()
            ax1.set_zorder(ax2.get_zorder() + 1)
            ax1.patch.set_visible(False)

            def draw(frame_end):
                ax1.clear()
                ax2.clear()
                end = int(frame_end)

                ax1.plot(time[:end], I1_pct[:end], label="Infected (Community 1)", color="C0", linestyle="--")
                ax1.plot(time[:end], I2_pct[:end], label="Infected (Community 2)", color="C1", linestyle="--")

                ax2.fill_between(time[:end], 0, C1_pct[:end], color="C0", alpha=0.1)
                ax2.fill_between(time[:end], 0, C2_pct[:end], color="C1", alpha=0.1)
                ax2.plot(time[:end], C1_pct[:end], color="C0", linestyle="-", linewidth=0.1)
                ax2.plot(time[:end], C2_pct[:end], color="C1", linestyle="-", linewidth=0.1)

                ax1.set_title(record["title"])
                ax1.set_xlabel(r"Scaled time ($\ln n$)")
                ax1.set_ylabel(r"Percentage of population (\%)")
                ax1.yaxis.set_label_position("left")
                ax1.yaxis.tick_left()
                ax1.set_xlim(0, animation_x_upper_lim)
                ax1.set_ylim(0, 30)
                ax2.set_ylabel(r"Cumulative percentage infected (\%)")
                ax2.yaxis.set_label_position("right")
                ax2.yaxis.tick_right()
                ax2.set_ylim(0, 30)
                ax1.legend(loc="upper right")
                fig.subplots_adjust(left=0.12, right=0.88, bottom=0.12, top=0.92)

            animation = FuncAnimation(fig, draw, frames=frame_indices, interval=1000 / ANIMATION_FPS, repeat=False)
            output_file = FIGURES_DIR / (
                f"simulation_n{record['n']}_b{record['beta']}_g{record['gamma']}_T{record['rho_T']}_H{record['rho_H']}_"
                f"c{record['c']}_t{record['max_time']}_seed{record['seed']}{record['suffix']}_linear_animated.gif"
            )
            animation.save(output_file, writer=PillowWriter(fps=ANIMATION_FPS))
            plt.close(fig)
            print(f"Saved {output_file}")

        from matplotlib.animation import FuncAnimation, PillowWriter

        animation_records = load_animation_records()
        animation_x_upper_lim = max(record["time"].max() for record in animation_records)
        y_upper_lim = 1e8

        for record in animation_records:
            time = record["time"]
            I1_log = np.where(record["I1"] > 0, record["I1"], np.nan)
            I2_log = np.where(record["I2"] > 0, record["I2"], np.nan)
            threshold_time = first_threshold_time(time, record["I1"] + record["R1"], record["n"]) if record["suffix"] else None
            frame_indices = np.unique(np.linspace(2, len(time), min(ANIMATION_FRAMES, len(time)), dtype=int))

            fig, ax = plt.subplots(figsize=(10, 6))

            def draw(frame_end):
                ax.clear()
                end = int(frame_end)
                current_time = time[end - 1]

                ax.plot(time[:end], I1_log[:end], label="Infected (Community 1)", color="C0", linestyle="-", linewidth=1.6)
                ax.plot(time[:end], I2_log[:end], label="Infected (Community 2)", color="C1", linestyle="-", linewidth=1.6)
                if threshold_time is not None and current_time >= threshold_time:
                    ax.axvline(threshold_time, color="red", linestyle="--", linewidth=1.0, label=r"$\tau_1(1\%)$")

                ax.set_title(record["title"])
                ax.set_yscale("log")
                ax.set_xlabel(r"Scaled time ($\ln n$)")
                ax.set_ylabel("Number of infectious individuals")
                ax.set_xlim(0, animation_x_upper_lim)
                ax.set_ylim(1, y_upper_lim * 1.2)
                ax.legend(loc="upper right")
                fig.subplots_adjust(left=0.14, right=0.96, bottom=0.12, top=0.92)

            animation = FuncAnimation(fig, draw, frames=frame_indices, interval=1000 / ANIMATION_FPS, repeat=False)
            output_file = FIGURES_DIR / (
                f"simulation_n{record['n']}_b{record['beta']}_g{record['gamma']}_T{record['rho_T']}_H{record['rho_H']}_"
                f"c{record['c']}_t{record['max_time']}_seed{record['seed']}{record['suffix']}_log_counts_animated.gif"
            )
            animation.save(output_file, writer=PillowWriter(fps=ANIMATION_FPS))
            plt.close(fig)
            print(f"Saved {output_file}")

if __name__ == "__main__":
    main()
