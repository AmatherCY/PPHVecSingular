# Analyzing Singular Patterns in Discrete Planar Vector Fields via Persistent Path Homology
![](https://github.com/AmatherCY/PPHVecSingular/blob/main/magfield.png)

This code is used to extract singularity locations from a discrete planar vector field using 1-dimensional persistent path homology and to analyze singular pattern changes in the time-varying vector field.

# Abstract
Analyzing singular patterns in vector fields is a fundamental problem in theoretical and practical domains due to the ability of such patterns to detect the intrinsic characteristics of vector fields. In this study, we propose an approach for analyzing singular patterns from discrete planar vector fields. Our method involves converting the planar discrete vector field into a specialized digraph and computing its one-dimensional persistent path homology. By analyzing the persistence diagram, we can determine the location of singularities, and the variations of singular patterns can also be analyzed. The experimental results demonstrate the effectiveness of our method in analyzing the singular patterns of noisy real-world vector fields and measuring the variations between different vector fields.

# Dependencies
This project requires the following Python packages:
- numpy
- pandas
- networkx
- gudhi

We recommend using Visual Studio Code to run this project.

## Environments
### Windows (Recommended)
Windows Users can install necessary Python packages using pip:
```bash
pip install numpy pandas networkx gudhi
```

### Linux
For Ubuntu 24.04 Users, due to Python environment protection ([PEP 668](https://peps.python.org/pep-0668/)), running `pip install` system-wide may result in an `externally-managed-environment` error. Please try the following method:

```bash
sudo apt update
sudo apt install python3-numpy python3-pandas python3-networkx python3-gudhi
```



# Data
The vector fields data is `.csv` format files, which is stored in `PPHVecSingular\IGRF\wfs`, `PPHVecSingular\2023kn\wfs`, `PPHVecSingular\2023sl\wfs`, and `PPHVecSingular\noise`.

# Detecting singularities
You can run the `Dip_pole_position.py` to track the geomagnetic poles example. The code will print the positions of geomagnetic poles in Fig.6 in the paper (or more specifically, the positioning data in Table.2)：
```bash
python Dip_pole_position.py
```

And you can run the `Typoon_trace.py` to track the typhoon centers example. By setting the variable `rig` as `2023kn` or `2023sl`, the code will print the centers of typhoon Kahnu (Fig.7 (a) in the paper) and Saola (Fig.7 (b) in the paper) respectively (or more specifically, the positioning data in Table.3 and 4).
```bash
python Typoon_trace.py
```

You can also run the `noise_compare.ipynb` to test the robustness using different noise level vector fields in `PPHVecSingular\noise` (Results of Fig.9 in the paper).

# Analyzing time-varying vector fields
After running `Typoon_trace.py` of two typhoon datasets, you can run the `VF_compare.py` to compute the bottleneck distance between persistence diagrams at neighboring moments (Results of Fig.8 in the paper), which can compare the changes of singular patterns of the wind fields. Similarly, you can set the variable `rig` as `2023kn` or `2023sl` for testing the two typhoons.
```bash
python VF_compare.py
```

# Acknowledgement
Our code for computing 1-d path homology is based on https://github.com/tianqicoding/1dPPH.
