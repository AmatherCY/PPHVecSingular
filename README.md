# Analyzing Singular Patterns in Discrete Planar Vector Fields via Persistent Path Homology
![](https://github.com/AmatherCY/PPHVecSingular/blob/main/magfield.png)

This code is used to extract singularity locations from a discrete planar vector field using 1-dimensional persistent path homology and to analyze singular pattern changes in the time-varying vector field.

# Dependencies
We use GUDHI to visualize persistence diagram:

``` python
pip install gudhi
```

# Abstract
Analyzing singular patterns in vector fields is a fundamental problem in theoretical and practical domains due to the ability of such patterns to detect the intrinsic characteristics of vector fields. In this study, we propose an approach for analyzing singular patterns from discrete planar vector fields. Our method involves converting the planar discrete vector field into a specialized digraph and computing its one-dimensional persistent path homology. By analyzing the persistence diagram, we can determine the location of singularities, and the variations of singular patterns can also be analyzed. The experimental results demonstrate the effectiveness of our method in analyzing the singular patterns of noisy real-world vector fields and measuring the variations between different vector fields.

# Data
The vector fields data is `.csv` format files, which is stored in `PPHVecSingular\IGRF\wfs`, `PPHVecSingular\2023kn\wfs`, `PPHVecSingular\2023sl\wfs`, and `PPHVecSingular\noise`.

# Detecting singularities
Please run the `Dip_pole_position.py` to track the geomagnetic poles example, and run the `Typoon_trace.py` to track the typhoon centers example.
You can also run the `noise_compare.ipynb` to test the robustness using different noise level vector fields in `PPHVecSingular\noise`.

# Analyzing time-varying vector fields
After running `Typoon_trace.py` of two typhoon datasets, you can run the `VF_compare.py` to compute the bottleneck distance between persistence diagrams at neighboring moments, which can compare the changes of singular patterns of the wind fields.

# Acknowledgement
Our code for computing 1-d path homology is based on https://github.com/tianqicoding/1dPPH.
