import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# === Parameters to customize ===
label_fontsize = 28
tick_fontsize = 24

x_tick_interval = 0.5
y_tick_interval_height = 0.05
y_tick_interval_1st = 0.5
y_tick_interval_2nd = 0.5

frame_interval = 1
start_frame = 9286 #9270
end_frame = 9286 #9270

file_path = "1_20250814_bar_load_100kHz_f_01_625n7_u2_0_10000.txt"

all_times = []
all_radii = []
all_heights = []
all_first_derivatives = []
all_second_derivatives = []

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

                dz_dr = []
                for i in range(len(r_sorted)):
                    if i == 0:
                        dz = z_sorted[i+1] - z_sorted[i]
                        dr = r_sorted[i+1] - r_sorted[i]
                    elif i == len(r_sorted) - 1:
                        dz = z_sorted[i] - z_sorted[i-1]
                        dr = r_sorted[i] - r_sorted[i-1]
                    else:
                        dz = z_sorted[i+1] - z_sorted[i-1]
                        dr = r_sorted[i+1] - r_sorted[i-1]
                    slope = dz / dr if dr != 0 else 0
                    dz_dr.append(slope)

                d2z_dr2 = []
                for i in range(len(r_sorted)):
                    if i == 0:
                        delta = dz_dr[i+1] - dz_dr[i]
                        dr = r_sorted[i+1] - r_sorted[i]
                    elif i == len(r_sorted) - 1:
                        delta = dz_dr[i] - dz_dr[i-1]
                        dr = r_sorted[i] - r_sorted[i-1]
                    else:
                        delta = dz_dr[i+1] - dz_dr[i-1]
                        dr = r_sorted[i+1] - r_sorted[i-1]
                    curvature = delta / dr if dr != 0 else 0
                    d2z_dr2.append(curvature)

                all_times.append([current_time] * len(r_sorted))
                all_radii.append(r_sorted)
                all_heights.append(z_sorted)
                all_first_derivatives.append(dz_dr)
                all_second_derivatives.append(d2z_dr2)

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

    if (
        current_time is not None and current_r and current_z and
        start_frame <= frame_number <= end_frame and
        (frame_number - start_frame) % frame_interval == 0
    ):
        sorted_pairs = sorted(zip(current_r, current_z))
        r_sorted, z_sorted = zip(*sorted_pairs)

        dz_dr = []
        for i in range(len(r_sorted)):
            if i == 0:
                dz = z_sorted[i+1] - z_sorted[i]
                dr = r_sorted[i+1] - r_sorted[i]
            elif i == len(r_sorted) - 1:
                dz = z_sorted[i] - z_sorted[i-1]
                dr = r_sorted[i] - r_sorted[i-1]
            else:
                dz = z_sorted[i+1] - z_sorted[i-1]
                dr = r_sorted[i+1] - r_sorted[i-1]
            slope = dz / dr if dr != 0 else 0
            dz_dr.append(slope)

        d2z_dr2 = []
        for i in range(len(r_sorted)):
            if i == 0:
                delta = dz_dr[i+1] - dz_dr[i]
                dr = r_sorted[i+1] - r_sorted[i]
            elif i == len(r_sorted) - 1:
                delta = dz_dr[i] - dz_dr[i-1]
                dr = r_sorted[i] - r_sorted[i-1]
            else:
                delta = dz_dr[i+1] - dz_dr[i-1]
                dr = r_sorted[i+1] - r_sorted[i-1]
            curvature = delta / dr if dr != 0 else 0
            d2z_dr2.append(curvature)

        all_times.append([current_time] * len(r_sorted))
        all_radii.append(r_sorted)
        all_heights.append(z_sorted)
        all_first_derivatives.append(dz_dr)
        all_second_derivatives.append(d2z_dr2)

# === Plot 1: Height (×0.1) ===
plt.figure(figsize=(10, 6))
for r, z in zip(all_radii, all_heights):
    r_shifted = [ri + 1 for ri in r]
    z_norm = [0.2 * zi for zi in z]
    plt.plot(r_shifted, z_norm, linewidth=5.0, color='navy')
plt.xlabel('Radial Distance (µm)', fontsize=label_fontsize)
plt.ylabel('Height (µm)', fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xlim(-1, 1)
plt.ylim(-0.01, 0.11)
plt.gca().xaxis.set_major_locator(MultipleLocator(x_tick_interval))
plt.gca().yaxis.set_major_locator(MultipleLocator(y_tick_interval_height))
plt.grid(True)
plt.tight_layout()
plt.show()

# === Plot 2: First Derivative (×0.1) ===
plt.figure(figsize=(10, 6))
for r, dz_dr in zip(all_radii, all_first_derivatives):
    r_shifted = [ri + 1 for ri in r]
    dz_dr_norm = [1 * val for val in dz_dr]
    plt.plot(r_shifted, dz_dr_norm, linewidth=5.0, color='navy')
plt.xlabel('Radial Distance (µm)', fontsize=label_fontsize)
plt.ylabel('First Derivative', fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xlim(-1, 1)
plt.ylim(-1.1, 1.1)
plt.gca().xaxis.set_major_locator(MultipleLocator(x_tick_interval))
plt.gca().yaxis.set_major_locator(MultipleLocator(y_tick_interval_1st))
plt.grid(True)
plt.tight_layout()
plt.show()

# === Plot 3: Second Derivative (×0.1) ===
plt.figure(figsize=(10, 6))
for r, d2z_dr2 in zip(all_radii, all_second_derivatives):
    r_filtered = []
    d2z_dr2_filtered = []
    for ri, val in zip(r, d2z_dr2):
        if not (-1.0006 <= ri <= -0.9994):  # filtering original r
            r_filtered.append(ri + 1)
            d2z_dr2_filtered.append(0.1 * val)
    plt.plot(r_filtered, d2z_dr2_filtered, linewidth=5.0, color='navy')
plt.xlabel('Radial Distance (µm)', fontsize=label_fontsize)
plt.ylabel('Second Derivative', fontsize=label_fontsize)
plt.xticks(fontsize=tick_fontsize)
plt.yticks(fontsize=tick_fontsize)
plt.xlim(-1, 1)
plt.ylim(-1.5, 1)
plt.gca().xaxis.set_major_locator(MultipleLocator(x_tick_interval))
plt.gca().yaxis.set_major_locator(MultipleLocator(y_tick_interval_2nd))
plt.grid(True)
plt.tight_layout()
plt.show()
