import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

# === Parameters ===
frame_interval = 1
start_frame = 15000
end_frame = 16000
total_frames = (end_frame - start_frame) // frame_interval + 1

file_path = "1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000.txt"
output_file = "1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000_2nd_derivative.csv"

# === Radial range and center
r_min_total = -2.0
r_max_total = 0.0
r_center = (r_min_total + r_max_total) / 2.0
cutoff_min = 0.0008  # µm
cutoff_max = 0.5   # µm

results = []
all_relative_positions = []
filtered_relative_positions = []

def compute_second_derivative(r_sorted, z_sorted):
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

    return d2z_dr2

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
                d2z_dr2 = compute_second_derivative(r_sorted, z_sorted)

                # --- Original max position
                abs_values = [abs(val) for val in d2z_dr2]
                max_index = abs_values.index(max(abs_values))
                r_max_position = r_sorted[max_index]
                max_magnitude = d2z_dr2[max_index]
                relative_position = abs(r_max_position - r_center)
                all_relative_positions.append(relative_position)

                # --- Filtered max position
                rel_pos_filtered = ''
                filtered_indices = [
                    i for i, r in enumerate(r_sorted)
                    if cutoff_min < abs(r - r_center) <= cutoff_max
                ]
                if filtered_indices:
                    filtered_values = [abs(d2z_dr2[i]) for i in filtered_indices]
                    max_i_filtered = filtered_indices[filtered_values.index(max(filtered_values))]
                    r_filtered = r_sorted[max_i_filtered]
                    rel_pos_filtered = abs(r_filtered - r_center)
                    filtered_relative_positions.append(rel_pos_filtered)

                results.append([current_time, max_magnitude, relative_position, rel_pos_filtered])

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

    # Final frame
    if (
        current_time is not None and current_r and current_z and
        start_frame <= frame_number <= end_frame and
        (frame_number - start_frame) % frame_interval == 0
    ):
        sorted_pairs = sorted(zip(current_r, current_z))
        r_sorted, z_sorted = zip(*sorted_pairs)
        d2z_dr2 = compute_second_derivative(r_sorted, z_sorted)

        abs_values = [abs(val) for val in d2z_dr2]
        max_index = abs_values.index(max(abs_values))
        r_max_position = r_sorted[max_index]
        max_magnitude = d2z_dr2[max_index]
        relative_position = abs(r_max_position - r_center)
        all_relative_positions.append(relative_position)

        rel_pos_filtered = ''
        filtered_indices = [
            i for i, r in enumerate(r_sorted)
            if cutoff_min < abs(r - r_center) <= cutoff_max
        ]
        if filtered_indices:
            filtered_values = [abs(d2z_dr2[i]) for i in filtered_indices]
            max_i_filtered = filtered_indices[filtered_values.index(max(filtered_values))]
            r_filtered = r_sorted[max_i_filtered]
            rel_pos_filtered = abs(r_filtered - r_center)
            filtered_relative_positions.append(rel_pos_filtered)

        results.append([current_time, max_magnitude, relative_position, rel_pos_filtered])

# === Save to CSV ===
with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow([
        'Time (s)',
        'Max 2nd Derivative (nm/um^2)',
        'Relative Position to Center (All Data) (um)',
        'Relative Position to Center (0.01 < |r-r₀| ≤ 0.5 µm) (um)'
    ])
    writer.writerows(results)

# === Plot Histograms ===
bin_size = 0.05
max_relative_pos = cutoff_max  # Only plot up to cutoff_max for filtered
bins = int(max_relative_pos / bin_size)


# --- Filtered histogram
counts_filt, bin_edges_filt = plt.hist(
    filtered_relative_positions,
    bins=bins,
    edgecolor='black'
)[0:2]

plt.clf()

normalized_counts_filt = counts_filt / total_frames
bin_centers_filt = 0.5 * (np.array(bin_edges_filt[:-1]) + np.array(bin_edges_filt[1:]))

plt.figure(figsize=(10, 6))
plt.bar(bin_centers_filt, normalized_counts_filt, width=0.9 * bin_size, edgecolor='black', color='blue')
plt.xlabel('Distance to Center (µm)', fontsize=28)
plt.ylabel('Fraction', fontsize=28)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

ax = plt.gca()
ax.text(0.98, 0.98, '1 kHz', transform=ax.transAxes,
        fontsize=28, verticalalignment='top', horizontalalignment='right')

plt.ylim(0, 0.7)
plt.xlim(0, 0.5)
plt.gca().xaxis.set_major_locator(MultipleLocator(0.1))
plt.gca().yaxis.set_major_locator(MultipleLocator(0.2))

plt.grid(True)
plt.tight_layout()
plt.savefig("1_2nd_hist_1kHz_2.png", dpi=300, bbox_inches=None, facecolor="white")
plt.show()