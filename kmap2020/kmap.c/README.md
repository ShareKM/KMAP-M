# Kinetic Modeling and Analysis Package (KMAP)

## Overview of the KMAP Library

The KMAP (Kinetic Modeling and Analysis Package) library provides a suite of code for implementing and applying various kinetic models to analyze dynamic positron emission tomography (PET) data. This library includes C/C++ code and MEX files that facilitate seamless integration with MATLAB, enabling efficient kinetic analysis in a flexible and user-friendly environment.

### Directory Structure

1. **`KMAP_source_files`**: Contains all the source files used for kinetic modeling, optimization, and common functions.
2. **`KMAP_mex_files`**: Contains the MEX files and their precompiled binaries for executing in MATLAB, along with all related code.

## Supported Kinetic Models

1. **One-Tissue Compartmental Model (1TCM)**
   - **Description**: A compartmental model where the tracer exchanges between plasma and a single tissue compartment.
   - **Parameters**: \( K_1 \), \( k_2 \), \( V_b \).

2. **Two-Tissue Compartmental Model (2TCM)**
   - **Description**: A more detailed model with two compartments representing different tissue types or states, providing a richer description of tracer kinetics.
   - **Parameters**: \( K_1 \), \( k_2 \), \( k_3 \), \( k_4 \), \( V_b \).
   - If \( k_4=0 \), this is an irreversible 2T model.

3. **Simplified Reference Tissue Model (SRTM)**
   - **Description**: This model is used when a reference region devoid of the receptor of interest is available, helping to estimate binding potential without requiring a blood input function.
   - **Parameters**: \( R_1 \), \( k_2 \), \( BP_{nd} \), \( V_b \).

4. **Liver Dual-Blood Input Function Model**
   - **Description**: A specialized model accounting for the dual-blood input function in liver tissue, where both hepatic artery and portal vein are considered.
   - **Parameters**: \( K_1 \), \( k_2 \), \( k_3 \), \( k_4 \), \( K_a \), \( f_a \), \( V_b \).

### For Further Instructions

To learn more about the detailed compilation instructions and how to use the files in this package, refer to the `README` inside the `KMAP_source_files` directory and `README` inside the `KMAP_mex_files` directory.

These files provide comprehensive guidance on compiling and executing the source and MEX files within MATLAB, along with additional details about the optimization methods employed.
