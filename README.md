# LiSBOA_scan_design
Tool for lidar scan design

# Description
This package allows the user to identify the optimal scan parameters for a lidar scan of different geometry. The applications are restricted to 2D (PPI or RHI) and volumetric scana that use azimuth and elevation with a Cartesian configuration. The limits of azimuth and elevation angles must also be symmetrical. Future releases may address more complex scan geometry.

The user defines the inputs in an Excel table. Then the LiSBOA_Pareto_design is run to calculate the three objective functions. Finally, the Pareto front and the full geometry of a selected scan can be visualized through the LiSBOA_design_postpro.

# Inputs
The inputs of the scan are defined in an Excel table like the data/Scan_info.xlsx. The full list of inputs is:
-Run [boolean]: if TRUE, the entry is processed, of FALSE, it is sckipped
- Scan name [string]: the name of the scan
- vec_dazi [deg]: the list of azimuth resolution to be tested
- vec_azi_max [deg]: the list of limits of the azimuth. The tested azimuth will go from -azi_max to +azi_max with a dazi step
- vec_dele [deg]: the list of elevation resolution to be tested
- vec_ele_max [deg]: the list of limits of the elevatin. The tested elevation will go from -ele_max to +ele_max with a dazi step
