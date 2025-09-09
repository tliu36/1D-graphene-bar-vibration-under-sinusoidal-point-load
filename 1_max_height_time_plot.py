import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import csv
from scipy.fft import fft, fftfreq

# === Parameters to customize ===
label_fontsize = 20
tick_fontsize = 16

frame_interval = 1      # Plot every Nth frame
start_frame = 0         # Minimum frame to include
end_frame = 10000        # Maximum frame to include
frame_rate = 100*5        # Frames per second (1 frame = 0.01s, so 100 Hz)
dt = 1 / frame_rate     # Time step

# === Load and parse the file ===
file_path = "1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000.txt"

all_times = []
all_radii = []
all_heights = []

with open(file_path, 'r') as f:
    current_time = None
    current_r = []
    current_z = []
    frame_number = None

    for line in f:
        parts = line.strip().split()

        if parts[0] == 'Frame':
            if (
                current_time is not None and current_r and current_z and
                start_frame <= frame_number <= end_frame and
                (frame_number - start_frame) % frame_interval == 0
            ):
                sorted_pairs = sorted(zip(current_r, current_z))
                r_sorted, z_sorted = zip(*sorted_pairs)

                all_times.append([current_time] * len(r_sorted))
                all_radii.append(r_sorted)
                all_heights.append(z_sorted)

            frame_number = int(float(parts[1]))
            current_time = frame_number / 100.0
            current_r = []
            current_z = []

        elif len(parts) == 3:
            try:
                r = float(parts[1])
                z = float(parts[2])
                current_r.append(r)
                current_z.append(z)
            except ValueError:
                continue

    # Final frame append
    if (
        current_time is not None and current_r and current_z and
        start_frame <= frame_number <= end_frame and
        (frame_number - start_frame) % frame_interval == 0
    ):
        sorted_pairs = sorted(zip(current_r, current_z))
        r_sorted, z_sorted = zip(*sorted_pairs)

        all_times.append([current_time] * len(r_sorted))
        all_radii.append(r_sorted)
        all_heights.append(z_sorted)

# === 2D Plot of Max Height vs. Time ===
frame_times = [t[0] for t in all_times]
max_heights = [max(z) for z in all_heights]

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(frame_times, max_heights, marker='o', linestyle='-', color='tab:blue')
ax.set_xlabel('Time (s)', fontsize=label_fontsize)
ax.set_ylabel('Max Height (µm)', fontsize=label_fontsize)
ax.tick_params(axis='x', labelsize=tick_fontsize)
ax.tick_params(axis='y', labelsize=tick_fontsize)
ax.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()

# === Save max height data to CSV ===
output_csv = "1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000_max_height_vs_time.csv"
with open(output_csv, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Time (s)', 'Max Height (µm)'])
    for t, h in zip(frame_times, max_heights):
        writer.writerow([t, h])
print(f"Saved max height vs. time data to: {output_csv}")

# === FFT Analysis of Max Height ===
N = len(max_heights)
yf = fft(max_heights)
xf = fftfreq(N, dt)[:N//2]

# === Plot FFT Spectrum ===
fig, ax_fft = plt.subplots(figsize=(10, 6))
ax_fft.plot(xf, 2.0 / N * np.abs(yf[0:N//2]), color='darkred')
ax_fft.set_xlabel('Frequency (Hz)', fontsize=label_fontsize)
ax_fft.set_ylabel('Amplitude', fontsize=label_fontsize)
ax_fft.tick_params(axis='x', labelsize=tick_fontsize)
ax_fft.tick_params(axis='y', labelsize=tick_fontsize)
ax_fft.grid(True, linestyle='--', alpha=0.6)
ax_fft.set_title('FFT of Max Height vs Time', fontsize=label_fontsize)
plt.tight_layout()
plt.show()
