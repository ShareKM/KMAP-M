# Kinetic Modeling and Analysis Package (KMAP)

## Overview

**KMAP** offers a collection of source code and wrapper functions designed to implement and apply different tracer kinetic models for analyzing dynamic positron emission tomography (PET) data, particularly in response to the challenges emerging in total-body PET kinetic modeling. The primary objective of this open-source package is to share tracer kinetic modeling techniques and offer kinetic modeling developers a foundation to build upon without starting from scratch. The package is originally developed at the University of California, Davis. Its open-source version is initiated as a part of the [Open Kinetic Modeling Initiative (OpenKMI)](https://www.openkmi.org/).

This package includes 
- `KMAP-C`: The C/C++ source code for implementing different kinetic modeling methods and the associated MEX files that facilitate seamless integration with MATLAB.
- `KMAP-M`: A MATLAB toolbox built on the KMAP-C for tracer kinetic modeling and can be independently used in MATLAB under Windows, Linux, and Macintosh systems. Recompiling may be needed if the pre-compiled files do not work.

## Ongoing Effort

This open-source package is provided "as is", without warranty. Our team is still iterating many functions in this package. 

## Licensing

This repository contains code that is licensed under [MIT License](KMAP-C/LICENSE).
