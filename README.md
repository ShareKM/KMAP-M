# Davis Kinetic Modeling and Analysis Package (DavisKMAP)

## Overview of DavisKMAP

DavisKMAP offers a collection of source code and wrapper functions designed to implement and apply different tracer kinetic models for analyzing dynamic positron emission tomography (PET) data, particularly considering the challenges emerging in total-body PET kinetic modeling. The primary objective of this open-source package is to share tracer kinetic modeling techniques and offer kinetic modeling developers a foundation to build upon without starting from scratch.

This package includes 
- `DavisKMAP-C`: The C/C++ source code for implementing different kinetic modeling methods and the associated MEX files that facilitate seamless integration with MATLAB.
- `DavisKMAP-M`: A MATLAB toolbox built on the DavisKMAP-C for tracer kinetic modeling and can be independently used in MATLAB under Windows, Linux, and Macintosh systems. Recompiling may be needed if the pre-compiled files do not work.

## Ongoing Effort

This open-source package is provided "as is", without warranty. Our team is still iterating many functions in this package. 

## Licensing

This repository contains code that is licensed under two different licenses:

- Code in `DavisKMAP-M` is licensed under [MIT License](DavisKMAP-M/LICENSE).
- Code in `DavisKMAP-C` is licensed under [Apache 2.0 License](DavisKMAP-C/LICENSE).
