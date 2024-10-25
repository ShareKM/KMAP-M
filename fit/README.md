# KMAP-M Fit Functions

This folder contains both original MATLAB code for kinetic fitting (e.g., the Patlak plot) and wrapper functions that are designed to interact with the compiled MEX binaries for kinetic modeling with compartmental modeling (e.g. the two-tissue model). These wrapper functions abstract the complexity of directly working with the MEX binaries, providing a more user-friendly interface for inputting data and obtaining results.

## Files in This Folder

- **kfit_1tcm.m**  
   Wrapper function for fitting time-activity curves (TACs) using the one-tissue compartment model (1TCM). This function calls the compiled MEX binary `kfit_1tcm_mex_omp` to perform the actual optimization.
   
- **kfit_2tcm.m**  
   Wrapper function for fitting TACs using the two-tissue compartment model (2TCM). This function calls the compiled MEX binary `kfit_2tcm_mex_omp` to perform the actual optimization.  

- **kfit_liver.m**  
   Wrapper function for fitting TACs using a liver model with seven parameters. This function calls the compiled MEX binary `kfit_liver_mex_omp` for optimization.

- **patlak_plot.m**
   A Matlab function that implements the Patlak Plot for estimating the influx rate parameter Ki and the intercept b for an irreversible tracer. The Patlak Plot requires a full-time input function.

- **relative_patlak.m**
   A Matlab function that implements the Replative Patlak plot method which does not require the early-time input function but only late time points.
  
- **run_ParametricImage.m**  
   This function is an example of how to generate parametric images with the kinetic models from the KMAP package. This function currently uses 'kfit_2tcm' to provide 0T/1T/2Ti models for TAC fitting but other models can also be included. 

   `Note`: This function requires an additional package 'PLOT_v1.0' which can be downloaded from https://wanglab.faculty.ucdavis.edu/code.
  The image data for testing this function can be obtained from https://drive.google.com/drive/folders/1OWG0FMlhuZL3vnv_fswXSaqXo8_762_5?usp=drive_link. 

## Usage Instructions

For example, to use the `kfit_2tcm.m` wrapper, you can call the function like this:

```matlab
[k, cfit] = kfit_2tcm(ct, cp, scant, k0, opt);
```
Specific examples of using many of the functions can be found in the `KMAP-M/demo` folder.

Ensure that the inputs are prepared according to the function specifications (refer to each function’s help documentation in MATLAB).

Before using any of those wrapper functions, ensure that the compiled MEX binaries (e.g., `kfit_2tcm_mex_omp`, `kfit_liver_mex_omp`) are available in your MATLAB path. These binaries should have been compiled using the appropriate source files in the `KMAP-C/mex` directory. If you wish to recompile the MEX binaries, refer to the `README.md` file inside the `KMAP-C/mex` folder for compilation instructions.

### Optional Settings (`opt`)

The `opt` structure used in those appropriate wrapper functions contains optional settings such as:

- `Decay`: Decay constant to be used in the model.
- `LowerBound`: Lower bounds for the parameters during optimization.
- `UpperBound`: Upper bounds for the parameters during optimization.
- `PrmSens`: A vector indicating which parameters are active during fitting.
- `MaxIter`: Maximum number of iterations for optimization.
- `TimeStep`: Time step to be used for numerical integration.

You can modify these settings by passing an `opt` structure when calling the wrapper functions. If no `opt` structure is provided, default values will be used.

### Error Handling

The wrapper functions check for NaN values, size mismatches between input data, and incorrect dimensions for scan times. If any of these issues are encountered, an error message will be displayed, and the fitting or TAC generation process will terminate.
