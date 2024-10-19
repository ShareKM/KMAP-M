# Open KMAP-M Toolbox

## Overview of OpenKMAP

The **Open Kinetic Modeling and Analysis Package (KMAP)** is an open-source software environment designed to implement and apply various tracer kinetic models for analyzing dynamic positron emission tomography (PET) data. It particularly focuses on addressing the challenges associated with total-body PET kinetic modeling, for example, considering time delay correction. The main goal of this open-source toolkit is to provide developers of tracer kinetic modeling with a foundational library to build upon, saving them from starting from scratch. The initial version of **KMAP** was developed at the University of California, Davis, and its open-source version was launched to support the [Open Kinetic Modeling Initiative (OpenKMI)](https://www.openkmi.org/).

## OpenKMAP-M Toolbox

**KMAP-M** is a MATLAB toolbox built on the C/C++ library of OpenKMAP for tracer kinetic modeling and can be independently used in MATLAB under Windows, Linux, and Macintosh systems. The toolbox also includes other kinetic modeling approaches that can be easily and efficiently implemented in MATLAB, such as the Patlak plot. 

The toolbox includes several folders:
- `demo`: Test and demonstration of how to use the fitting code

- `fit`: Containing MATLAB code and wrapper functions for kinetic modeling

- `tac`: Containing code for generating time activity curves

- `data`: Data used for testing and demonstrating the code in this toolbox

- `utils`: Utility functions 

- `Precompiled_Binaries`: precompiled binaries for TAC fitting with compartmental models. Recompiling using the [KMAP-C](https://github.com/openkmi/KMAP-C) Toolkit may be needed if the precompiled files do not work for your system.

## Ongoing Effort

This open-source package is provided "as is", without warranty. Our team is still iterating many functions in this package. 

## Licensing

This repository contains code that is licensed under [MIT License](KMAP-C/LICENSE).
