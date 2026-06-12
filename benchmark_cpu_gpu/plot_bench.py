import os
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Configuration ─────────────────────────────────────────────────────────────

# Backends to plot (one subplot each); edit order freely
BACKENDS = ["cpu", "gpu", "cuda", "hip"]

MACHINES = ["local","sandbox"] # "bigfoot_NV100","JZ_NV100"
GPU_NAMES = {"local":"RTX1000",
             "sandbox":"RTX4090",
             "bigfoot_V100":"V100"}
CPU_NAMES = {"local":"i7x8",
             "sandbox":"i7x8"}

# Colour palette (one colour per machine, cycles if more machines than colours)
PALETTE = plt.cm.tab10.colors

BASELINE_MACHINE = "sandbox"
BASELINE_BACKEND = "gpu" # == openGL

FILES = {'local':'results_bench_local_turbulence',
         'sandbox':'results_bench_sandbox_turbulence'}

CASE = 'turbulence'

# ── Parse files ───────────────────────────────────────────────────────────────

def check_line_speed(line: str) -> float | None:
    if line.startswith("#"):
        m = re.search(r"([\d.eE+\-]+)\s+points\.step/s", line)
        if m:
            return float(m.group(1))
        else:
            return None


def parse_file(filepath: str) -> dict:
    results = {}

    with open(filepath, 'r') as f:
        lines = f.readlines()

    current_backend = None
    current_resolution = None

    for line in lines:
        line = line.strip()
        
        # Match backend (cpu/gpu)
        if line in ('cpu', 'gpu','cuda','hip'):
            current_backend = line
            results[current_backend] = {}
        
        # Match resolution (standalone integer)
        elif re.fullmatch(r'\d+', line):
            current_resolution = int(line)
        
        # Match the summary comment line with speed
        elif line.startswith('#') and 'points.step/s' in line:
            match = re.search(r'([\d.e+]+)\s+points\.step/s', line)
            if match and current_backend and current_resolution is not None:
                results[current_backend][current_resolution] = float(match.group(1))

    return results

#print(parse_file("./results_bench_local_turbulence"))

data = {}
for machine in MACHINES:
    data[machine] = parse_file(FILES[machine])

fig, axes = plt.subplots(2, len(BACKENDS), figsize=(9, 6), constrained_layout=True)
machine_colors = {m: PALETTE[i % len(PALETTE)] for i, m in enumerate(MACHINES)}
baseline = data[BASELINE_MACHINE][BASELINE_BACKEND]
baseline_speeds = list(baseline.values())

all_res = [float(value) for value in baseline.keys()]
n_m   = len(data.keys())
x     = np.arange(len(all_res))
width = 0.8 / n_m

maxspeed = 0.
for i, machine in enumerate(data.keys()):
    for j, backend in enumerate(data[machine].keys()):
        ax_speed = axes[0, j]
        ax_speedup = axes[1, j]
        # ── speed ──
        speeds = list(data[machine][backend].values())
        offset  = (i - n_m / 2 + 0.5) * width
        x = np.arange(len(data[machine][backend].keys()))
        if backend == "cpu":
            label=CPU_NAMES[machine]
        else:
            label=GPU_NAMES[machine]
        bars    = ax_speed.bar(
            x + offset, speeds, width,
            label=label,
            color=machine_colors[machine], edgecolor="white", linewidth=0.5,
        )
        ax_speed.set_title(backend)
        for speed in speeds:
            if speed>maxspeed:
                maxspeed = speed

        # ── speedup ──
        speedups = np.zeros(len(speeds))
        for k in range(len(speeds)):
            speedups[k] = speeds[k]/baseline_speeds[k]
        offset  = (i - n_m / 2 + 0.5) * width
        x = np.arange(len(data[machine][backend].keys()))
        bars    = ax_speedup.bar(
            x + offset, speedups, width,
            label=machine,
            color=machine_colors[machine], edgecolor="white", linewidth=0.5,
        )

for backend in range(len(BACKENDS)):
    axes[0,backend].set_ylim([0, 1.1*maxspeed])
    axes[1,backend].set_ylim([0,1.1])
axes[0,0].set_ylabel('speed (points.step/s)')
axes[0,0].legend(fontsize=8, framealpha=0.8)
axes[0,1].legend(fontsize=8, framealpha=0.8)
axes[1,0].set_ylabel('speedup (vs RTX4090)')

fig.savefig(f"bench_{CASE}.png", dpi=150) #, bbox_inches="tight")
plt.show()


