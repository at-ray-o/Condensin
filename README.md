# Condensin Simulation Package

Contains MD-simulation code written to simulate our model of Condensin using OpenMM.

## Description

This python script can be used to simulate equilibirum properties of Condensin and DNA by changing the number of steps Condensin spends in either B-state or O-state. The MD simulation is written using the OpenMM package. The script takes bond-length, bond-angle and initial particle positions as input. Finally it produces a trajectory and topology file as an output. The trajectory is written in dcd format and can be read using VMD.

## Documentation

## Requirements
 - Basic underderstanding of bash, python and how to use terminal
 - Python 3.7+
 - Python package installer conda, to install OpenMM. 
 - Linux, Mac or Windown OS

## Installing OpenMM
 - For general OpenMM installation instructions see http://docs.openmm.org/latest/userguide/application.html#installing-openmm. 
 - We used `conda install -c conda-forge openmm cudatoolkit==10.0`
 - Clone/download the repository
 - Move to the downloaded folder.
 - Run the following script in a bash window capable of running python, `python CG_model_Condensin_Elbow.py `

## Input Files
 - `coord_start.pdb` - Contains initial starting structure
 - `bondList.dat` - Contains list of bonds and bond-lengths
 - `DNAangles.dat` - Contains list of bond-angles used in DNA
 - `CondensinAngles.dat` - Contains list of angles used in Condensin model
 - `particleList.dat` - Contains particle names and sizes

## Output Files
 - `topology.pdb` - Contains a topology file to be used with dcd file
 - `output.dcd` - Contains output trajectory in dcd format
 - The bash window will also print all the energy values for each frame in the trajectory. This output can be redirected to a file if necessary.

## Examples

Example Input:
`python CG_model_Condensin_Elbow.py`

Example Output:

`Loaded everything
Setup forces

#"Step"   "Time (ps)"   "Potential Energy (kJ/mole)"   "Total Energy (kJ/mole)"
10000   39.99999999999959   710.12109375   936.8952330127358
20000   80.00000000000584   754.1318359375   1012.8524764608592
30000   120.00000000005473   790.8985595703125   1044.234991952777
40000   159.9999999999899   748.962646484375   1043.2672696188092
50000   199.9999999998967   765.4717407226562   1036.2402911083773
60000   239.99999999980346   747.9782104492188   1005.5008238162845
70000   279.99999999988074   767.0703125   1018.1886996394023
80000   320.00000000007174   796.067138671875   1042.7340556159616
90000   360.00000000026273   776.9962158203125   1050.8454275801778
100000   400.0000000004537   768.4482421875   1073.195808040211
110000   440.0000000006447   1028.3388671875   1390.8749343901873
120000   480.0000000008357   955.2681274414062   1287.7714983522892
130000   520.0000000010267   959.2990112304688   1271.9728351458907
140000   560.0000000012177   949.7865600585938   1213.0935844928026
150000   600.0000000014087   916.6130981445312   1194.2609043851262
160000   640.0000000015997   926.5377197265625   1191.74046555534
170000   680.0000000017907   960.9150390625   1249.244273437187
180000   720.0000000019817   947.8404541015625   1243.5862293243408
190000   760.0000000021727   926.4295654296875   1208.819168381393
200000   800.0000000023637   917.871826171875   1203.7580848727375
210000   840.0000000025547   811.3917236328125   1119.5756307952106
220000   880.0000000027457   782.7474365234375   1067.1339964009821
230000   920.0000000029366   787.230712890625   1054.1576135675423
240000   960.0000000031276   806.4061889648438   1092.6120233759284
250000   1000.0000000033186   789.201171875   1084.6818463131785
260000   1040.0000000030548   796.983642578125   1079.317869907245
270000   1080.000000002109   770.011474609375   1036.6852607652545
280000   1120.000000001163   801.9671630859375   1060.0420704837888
290000   1160.0000000002171   750.493896484375   1014.5487550422549
300000   1199.9999999992713   761.3817749023438   1019.675249459222
Time Elapsed :  41.842265129089355
`
