# PPHVecSingular

This code is used to extract singularity locations from a discrete planar vector field using 1-dimensional persistent path homology and to analyze singular pattern changes in the time-varying vector field.

# Data

The vector fields data is the csv. format files, which is stored in `PPHVecSingular\IGRF\wfs`, `PPHVecSingular\2023kn\wfs`, `PPHVecSingular\2023sl\wfs`, and `PPHVecSingular\noise`.

# Detecting singularities
Please run the `Dip_pole_position.py` to track the geomagnetic poles example, and run the `Typoon_trace.py` to track the typhoon centers example.
You can also run the `noise_compare.ipynb` to test the robustness using different noise level vector fields in `PPHVecSingular\noise`.

# Analyzing time-varying vector fields
After running `Typoon_trace.py` of two typhoon datasets, you can run the `VF_compare.py` to compute the bottleneck distance between persistence diagrams at neighboring moments, which can compare the changes of singular patterns of the wind fields.

# Acknowledgement
Our code for computing 1-d path homology is based on https://github.com/tianqicoding/1dPPH.
