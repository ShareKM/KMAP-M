# C-KMAP Toolkit

## Overview

**Open Kinetic Modeling and Analysis Package (OpenKMAP)** offers a developmental environment to implement and apply different tracer kinetic models for analyzing dynamic positron emission tomography (PET) data, particularly to address the challenges emerging in total-body PET kinetic modeling. The primary objective of this open-source package is to provide tracer kinetic modeling developers with a core library to build upon without starting from scratch. The early version of the **KMAP** package was originally developed at the University of California, Davis. Its open-source version is initiated as a part of the [Open Kinetic Modeling Initiative (OpenKMI)](https://www.openkmi.org/).

The **C-KMAP Toolkit** is a core C/C++ source code library to define and implement the kinetic models, optimization algorithms, and utility functions. For parametric imaging (voxel-wise kinetic modeling), it adopts the OpenMP parallel programming model for acceleration. The package also includes wapper functions designed to integrate the kinetic modeling implementations in MATLAB.

## Licensing

This repository is licensed under [MIT License](KMAP-C/LICENSE).
