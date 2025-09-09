This is an example ABAQUS simulation setup for a thin 1D graphene bar's vibration under sinusoidal point load in the center of the bar with two ends fixed.

1. The ABAQUS simulation setup for this iD model are in: 20250814_bar_load_1_10_100_200kHz_f.cae and 20250814_bar_load_1_10_100_200kHz_f.jnl   (The simulation results file is not included here as it is too large, but can be obtained by running the setup by choosing the desired frequency and set related time increment for that frequency.)

2. Once you have your OBD file from running the simulation setup, the displacement of each mesh of the bar can be processed by python code: 1_height_position_time_ave.py, all information are saved into: 1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000.txt (updated the name based on your obd file and the fine name of your uotput)

3. Then the max height of the center can be processed by the python code: 1_max_height_time_plot.py (This is used to check if steady state is achieved)

4. The evolution of the bar from the largest displacement to the minimum displacement can be obtained by the python code: 1_height_position_time_peak_frequency_3d_1kHz.py (Fig. 4)

5. The pattern, ist derivative, and 2nd derivative of slecyed frames can be obtained by python code: 1_height_position_1st_2nd_2d_100kHz.py

6: The histogram plot of the position of maximum 2nd derivative can be btained by python code: 1_height_2nd_max_position_new_1kHz_2.py

7. Example vibration pattern of 1 kHz and 200 kHz can be seen in 1_20250814_bar_load_1kHz.avi and 1_20250814_bar_load_200kHz