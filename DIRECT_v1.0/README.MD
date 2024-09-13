# KMAP Wrapper Functions for MEX Binaries

This folder contains MATLAB wrapper functions that are designed to interact with the compiled MEX binaries for kinetic modeling. These functions abstract the complexity of directly working with the MEX binaries, providing a more user-friendly interface for inputting data and obtaining results.

## Files in This Folder

1. **kfit_2tcm.m**  
   Wrapper function for fitting time-activity curves (TACs) using the two-tissue compartment model (2TCM). This function calls the compiled MEX binary `kfit_2tcm_mex_omp` to perform the actual optimization.  
   - **Inputs**:  
     - `ct`: Measured time-activity curve (TAC) data.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `k0`: Initial kinetic parameters.
     - `opt`: Fitting options structure.
     - `wt`: Frame weights (optional).
     - `cwb`: Whole blood data (optional).
   - **Outputs**:  
     - `k`: Fitted kinetic parameters.
     - `cfit`: Fitted TAC curve.

2. **kfit_liver.m**  
   Wrapper function for fitting TACs using a liver model with seven parameters. This function calls the compiled MEX binary `kfit_liver_mex_omp` for optimization.  
   - **Inputs**:  
     - `ct`: Measured TAC data.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `k0`: Initial kinetic parameters.
     - `opt`: Fitting options structure.
     - `wt`: Frame weights (optional).
     - `cwb`: Whole blood data (optional).
   - **Outputs**:  
     - `k`: Fitted kinetic parameters.
     - `cfit`: Fitted TAC curve.
     - `ctis`: Tissue response function (optional).

3. **ktac_2tcm.m**  
   Wrapper function for generating TACs using the two-tissue compartment model (2TCM). This function calls the compiled MEX binary `ktac_2tcm_mex` for TAC computation.  
   - **Inputs**:  
     - `k`: Kinetic parameters.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `opt`: Optional settings structure.
     - `cwb`: Whole blood concentration (optional).
   - **Outputs**:  
     - `ct`: Generated time-activity curve.
     - `st`: Sensitivity matrix (optional).

4. **ktac_liver.m**  
   Wrapper function for generating TACs from kinetic parameters for liver data. This function calls the compiled MEX binary `ktac_liver_mex` for TAC computation.  
   - **Inputs**:  
     - `k`: Kinetic parameters.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `opt`: Optional settings structure.
     - `cwb`: Whole blood concentration (optional).
   - **Outputs**:  
     - `ct`: Generated time-activity curve.
     - `st`: Sensitivity matrix (optional).

## Usage Instructions

Before using any of these wrapper functions, ensure that the compiled MEX binaries (e.g., `kfit_2tcm_mex_omp`, `kfit_liver_mex_omp`, `ktac_2tcm_mex`, and `ktac_liver_mex`) are available in your MATLAB path. These binaries should have been compiled using the appropriate source files in the `kmap2020/KMAP_mex_files` directory.

**Note:** The folder `kmap2020/KMAP_mex_files` contains both the MEX files and their precompiled binaries for different operating systems.  
- **Precompiled Binaries**: The precompiled MEX binaries for Windows, Linux, and macOS are organized in separate subdirectories within `kmap2020/KMAP_mex_files`. Use the binaries corresponding to your system architecture (`*.mexw64`, `*.mexa64`, or `*.mexmaci64`).
- If you wish to recompile the MEX binaries, refer to the `README.md` file inside the `KMAP_mex_files` folder for compilation instructions.

For example, to use the `kfit_2tcm.m` wrapper, you can call the function like this:

```matlab
[k, cfit] = kfit_2tcm(ct, cp, scant, k0, opt);
```

Ensure that the inputs are prepared according to the function specifications (refer to each function’s help documentation in MATLAB).

### Optional Settings (`opt`)

The `opt` structure used in these wrapper functions contains optional settings such as:

- `Decay`: Decay constant to be used in the model.
- `LowerBound`: Lower bounds for the parameters during optimization.
- `UpperBound`: Upper bounds for the parameters during optimization.
- `PrmSens`: A vector indicating which parameters are active during fitting.
- `MaxIter`: Maximum number of iterations for optimization.
- `TimeStep`: Time step to be used for numerical integration.

You can modify these settings by passing an `opt` structure when calling the wrapper functions. If no `opt` structure is provided, default values will be used.

### Error Handling

The wrapper functions check for NaN values, size mismatches between input data, and incorrect dimensions for scan times. If any of these issues are encountered, an error message will be displayed, and the fitting or TAC generation process will terminate.
