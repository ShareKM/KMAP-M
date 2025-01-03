<div align="center">
  <img src="https://github.com/user-attachments/assets/5079dd48-e136-466c-95ca-14155b15781d" width="400" >
</div>

## Overview of KMAP

The **Kinetic Modeling and Analysis Package (KMAP)** is an open-source software environment designed to implement and apply various tracer kinetic models for analyzing dynamic positron emission tomography (PET) data. It particularly focuses on addressing the challenges associated with total-body PET kinetic modeling, for example, considering time delay correction and other factors. The main goal of this open-source package is to provide developers of tracer kinetic modeling with a foundational library to build upon, saving them from starting from scratch. The initial version of **KMAP** was developed at the University of California, Davis, and its open-source version was launched to support the [Open Kinetic Modeling Initiative](https://www.openkmi.org/).

## The KMAP-M Toolbox

**KMAP-M** is a MATLAB toolbox built on the C/C++ library of **KMAP** for tracer kinetic modeling and can be independently used in MATLAB under Windows, Linux, and Macintosh systems. The toolbox also includes other kinetic modeling approaches that can be efficiently implemented in MATLAB (e.g., the Patlak plot). 

The toolbox includes several folders:
- `demo`: Test and demonstration of how to use the fitting code

- `fit`: Containing MATLAB code and wrapper functions for kinetic modeling

- `tac`: Containing code for generating time activity curves

- `data`: Data used for testing and demonstrating the code in this toolbox

- `utils`: Utility functions 

- `Precompiled_Binaries`: precompiled binaries for TAC fitting with compartmental models. Recompiling using the [**KMAP-C**](https://github.com/sharekm/KMAP-C) Toolkit may be needed if the precompiled files do not work for your system.

## MATLAB Path

Before using this toolbox, ensure that all the files are available in your MATLAB path. This can be done by running the path setup file `setup.m` under the main directory of KMAP-M.

>[!NOTE]
>**Ongoing Effort**: This toolbox is provided “as is” without warranty. The development team is continually updating the toolbox.

## Contributors

See the full list of [Contributors](CONTRIBUTORS.md).

## Licensing

This repository contains code that is licensed under [MIT License](KMAP-C/LICENSE).
