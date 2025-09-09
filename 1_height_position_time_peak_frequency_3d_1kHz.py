import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D

# === Parameters to customize ===
label_fontsize = 28
tick_fontsize = 24

x_tick_interval = 0.5  # time
y_tick_interval = 0.5  # radial distance (normalized)
z_tick_interval = 0.01  # height

frame_interval = 1
start_frame = 14004
end_frame = 14030
time_per_frame = 0.0625  # ms per frame
reference_frame = start_frame  # set this as time zero

# === View angle customization for 3D plot ===
elev_angle = 20  # elevation angle
azim_angle = -20  # azimuth angle

# === Radial range and normalization ===
r_min_total = -2
r_max_total = 0.0
focus_fraction = 0.3
midpoint = (r_max_total + r_min_total) / 2.0
range_half_width = (r_max_total - r_min_total) * focus_fraction / 2.0
focus_min = midpoint - range_half_width
focus_max = midpoint + range_half_width

# === Thresholds for noise
single_peak_slope_threshold = 0.000005
multi_peak_slope_threshold = 0.0000005

file_path = "1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000.txt"

all_times = []
all_radii = []
all_heights = []

single_peak_count = 0
multiple_peak_count = 0
total_frames = 0

def is_monotonic(seq, increasing=True, slope_thresh=0.0005):
    if len(seq) < 2:
        return False
    if increasing:
        return all((y - x) >= -slope_thresh for x, y in zip(seq, seq[1:]))
    else:
        return all((x - y) >= -slope_thresh for x, y in zip(seq, seq[1:]))

def clean_sign(slope, slope_thresh=0.0005):
    if abs(slope) < slope_thresh:
        return 0
    return 1 if slope > 0 else -1

def count_strict_reversals(signs):
    clean = [s for s in signs if s != 0]
    return sum(1 for a, b in zip(clean[:-1], clean[1:]) if a != b)

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

                normalized_time = (frame_number - reference_frame) * time_per_frame
                all_times.append([normalized_time] * len(r_sorted))
                all_radii.append(r_sorted)
                all_heights.append(z_sorted)
                total_frames += 1

                focus_data = [(r, z) for r, z in zip(r_sorted, z_sorted)
                              if focus_min <= r <= focus_max and z > 0]
                if len(focus_data) >= 3:
                    r_focus, z_focus = zip(*focus_data)
                    mid = (focus_min + focus_max) / 2

                    left = [z for r, z in focus_data if r <= mid]
                    right = [z for r, z in focus_data if r > mid]

                    if is_monotonic(left, increasing=True, slope_thresh=single_peak_slope_threshold) and \
                       is_monotonic(right, increasing=False, slope_thresh=single_peak_slope_threshold):
                        single_peak_count += 1

                    slopes = [z2 - z1 for z1, z2 in zip(z_focus[:-1], z_focus[1:])]
                    signs = [clean_sign(s, slope_thresh=multi_peak_slope_threshold) for s in slopes]
                    reversals = count_strict_reversals(signs)

                    if reversals >= 2:
                        multiple_peak_count += 1

            frame_number = int(float(parts[1]))
            current_time = frame_number * time_per_frame
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

        normalized_time = (frame_number - reference_frame) * time_per_frame
        all_times.append([normalized_time] * len(r_sorted))
        all_radii.append(r_sorted)
        all_heights.append(z_sorted)
        total_frames += 1

        focus_data = [(r, z) for r, z in zip(r_sorted, z_sorted)
                      if focus_min <= r <= focus_max and z > 0]
        if len(focus_data) >= 3:
            r_focus, z_focus = zip(*focus_data)
            mid = (focus_min + focus_max) / 2

            left = [z for r, z in focus_data if r <= mid]
            right = [z for r, z in focus_data if r > mid]

            if is_monotonic(left, increasing=True, slope_thresh=single_peak_slope_threshold) and \
               is_monotonic(right, increasing=False, slope_thresh=single_peak_slope_threshold):
                single_peak_count += 1

            slopes = [z2 - z1 for z1, z2 in zip(z_focus[:-1], z_focus[1:])]
            signs = [clean_sign(s, slope_thresh=multi_peak_slope_threshold) for s in slopes]
            reversals = count_strict_reversals(signs)

            if reversals >= 2:
                multiple_peak_count += 1

# === Output summary ===
print(f"Total frames processed: {total_frames}")
print(f"Single peak frames: {single_peak_count}")
print(f"Multiple peak (double/triple) frames: {multiple_peak_count}")

# === 3D Plotting (time = X, radial = Y, height = Z) ===
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

for time_series, r_list, z_list in zip(all_times, all_radii, all_heights):
    time_values = time_series  # Already normalized from 0 to 0.6 s
    r_normalized = [(2 * (r - r_min_total) / (r_max_total - r_min_total)) - 1 for r in r_list]
    z_scaled = [0.002 * z for z in z_list]  # Height normalized
    ax.plot(time_values, r_normalized, z_scaled, linewidth=4.0)

# Labels
ax.set_xlabel('Time (ms)', fontsize=label_fontsize, labelpad=20)
ax.set_ylabel('Radial Distance (µm)', fontsize=label_fontsize, labelpad=20)
ax.set_zlabel('Height (µm)', fontsize=label_fontsize, labelpad=20)

# Tick font size
ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

# Axis range
ax.set_xlim([0.0, 1.7])
ax.set_ylim([-1.1, 1.1])
ax.set_zlim([-0.11, 0.11])
ax.xaxis.set_major_locator(MultipleLocator(0.5))   # X ticks every 0.2
ax.yaxis.set_major_locator(MultipleLocator(0.5))   # Y ticks every 0.5
ax.zaxis.set_major_locator(MultipleLocator(0.05))  # Z ticks every 0.05

# === View perspective ===
ax.view_init(elev=elev_angle, azim=azim_angle)
#ax.text(0.2, -3, -0.1, "0.1 kHz", fontsize=28, color='black')


# Save with the intended canvas ratio (avoid 'tight' which can crop)
plt.savefig("1_displacement_1kHz.png", dpi=300, bbox_inches=None, facecolor="white")
plt.show()