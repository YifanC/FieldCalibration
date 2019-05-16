# FieldCalibration

To execute spatial displacment calculation with laser tracks

```$ ./FieldCal -D -d [number of track subsets] -N [number of iteration steps] -t -A [filename1] [filename2]... ```

// "-D": one must have it to run the spatial displacement calculation

// with "-C": correction map (reco -> true)

// without "-C": distortion map (true -> reco)

// "-t": by default use two side iteration

// "-A": by default use boundary condition


```$ ./FieldCal -E -T -M [number of toy throw] [filename] ```

// "-E": one must have it to run the E-field calculation 

// "-T": enable toy throw to access the error of E-field calculation

// "-M": number of toy throw (number of generated displacement map)

// Input TH3(s) from the root file are expected to be called "Reco_Displacement_X", "Reco_Displacement_Y", "Reco_Displacement_Z"


One may improve the way of command line execution experience for the users of the code.

To get the latest working version, please check out "develop" branch
