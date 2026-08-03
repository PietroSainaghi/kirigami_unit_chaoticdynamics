# kirigami_unit_chaoticdynamics

## Dataset S1

Data for simulation parametric resonance study, as text files. A CSV key file indicates the amplitude and frequencies used in files labelled with the corresponding row number (starting from zero). Files labelled as ‘Disp’ include the time steps for each dataset, files laballed ‘CenterLeft’ indicate the displacements in the probe node.

## Dataset S2

Data for experimental parametric resonance study, in Moku and MATLAB data formats.  Each subfolder indicates the sample used to collect the data. Each file name indicates the sample number, amplifier gain, function generator frequency, vibrometer gain for the corresponding data. The data in each file is stored in a table, where the first column indicates time, the second column indicates voltage reading from the impedance head, the third column indicates voltage reading from the vibrometer.

## Dataset S3

Data for key recognition study. Each file name indicates the sample identifier (A through E), amplifier gain, function generator frequency, vibrometer gain for the corresponding data. The data in each file is stored in a table, where the first column indicates time, the second column indicates voltage reading from the impedance head, the third column indicates voltage reading from the vibrometer.

## Software S1

Python scripts used to generate Abaqus finite-element simulation environment. Script 1 generates the first buckling mode displacement field, and computes the buckling load for a given geometry, determined by columns A through P in the ‘KirigamiFrequencyData_Damped’ CSV file. Script 2 computes the linear resonance eigenmodes for a given geometry. The results from 1 and 2 are used to fill in the values in columns R (buckling load), S (fraction of buckling load), T (load frequency) in the CSV. Script 3 runs dynamic simulations for a subset of the parametric resonance space determined by the selected rows of the CSV file.

## Software S2

MATLAB program for key recognition algorithm used in the feasibility study. The program follows the logic as highlighted in figure S5. The main script ‘Frugal_PhysicallyUniqueKey’ computes the metric vectors from each experimental test, and performs the testing and training study. The function ‘test_train_PhysicallyUnique’ contains the objective function computed every iteration of the study.

