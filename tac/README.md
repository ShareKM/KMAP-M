# KMAP-M Functions for Generating TACs

This folder contains MATLAB wrapper functions for generating time activity curves (TACs) and optional sensitivity functions from kinetic parameters.

## Files in This Folder

- **ktac_1tcm.m**  
   Wrapper function for generating TACs using the one-tissue compartment model (1TCM). This function calls the compiled MEX binary `ktac_1tcm_mex` for TAC computation.  
   - **Inputs**:  
     - `k`: Kinetic parameters.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `opt`: Optional settings structure.
     - `cwb`: Whole blood concentration (optional).
   - **Outputs**:  
     - `ct`: Generated time-activity curve.
     - `st`: Sensitivity functions (optional).

- **ktac_2tcm.m**  
   Wrapper function for generating TACs using the two-tissue compartment model (2TCM). This function calls the compiled MEX binary `ktac_2tcm_mex` for TAC computation.  
   - **Inputs**:  
     - `k`: Kinetic parameters.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `opt`: Optional settings structure.
     - `cwb`: Whole blood concentration (optional).
   - **Outputs**:  
     - `ct`: Generated time-activity curve.
     - `st`: Sensitivity functions (optional).

- **ktac_liver.m**  
   Wrapper function for generating TACs from kinetic parameters for liver data. This function calls the compiled MEX binary `ktac_liver_mex` for TAC computation.  
   - **Inputs**:  
     - `k`: Kinetic parameters.
     - `cp`: Plasma input function.
     - `scant`: Scan time data (start and end times in seconds).
     - `opt`: Optional settings structure.
     - `cwb`: Whole blood concentration (optional).
   - **Outputs**:  
     - `ct`: Generated time-activity curve.
     - `st`: Sensitivity functions (optional).
