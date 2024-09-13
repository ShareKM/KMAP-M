# KMAP MEX Files for MATLAB

This folder contains MEX files, along with their source code, that allow for running the KMAP library in MATLAB. These files enable fitting and generating time activity curves for kinetic models using the Levenberg-Marquardt algorithm and other advanced optimization methods. Precompiled MEX binaries are also included and can be used directly in MATLAB without recompilation, provided they are compatible with your system architecture.

Each MEX file is accompanied by detailed instructions on its usage, input, output, and compilation. Please refer to the comments within each MEX source file (`*.cpp`) for precise information on how to work with these files.

## Files in This Folder

- **Precompiled MEX Files**: Precompiled MEX binaries for Windows (`*.mexw64`), Linux (`*.mexa64`), and macOS (`*.mexmaci64`) are available in three separate folders: `Precompiled_Windows_Binaries`, `Precompiled_Linux_Binaries`, and `Precompiled_MacOS_Binaries`. Go to the respective folder depending on your system to access the binaries for direct use in MATLAB.
- **Source Code Files**: MEX source files (`*.cpp`) that implement the kinetic models in conjunction with the core source files from the `KMAP_source_files` folder.

### Compilation Instruction:

The compilation command should be run in the MATLAB Command Window or a terminal/shell where you have access to the MATLAB `mex` compiler. This command compiles the C++ source code into a MEX file that can be executed within MATLAB.

**Important:** Before compiling any of the MEX files, ensure that all required source files (`kinlib.h`, `kinlib_common.cpp`, `kinlib_models.cpp`, `kinlib_optimization.cpp`) from the `KMAP_source_files` directory are placed in the same directory as the MEX files.

### Usage:

The usage command should be run within the MATLAB environment after the MEX file has been successfully compiled. It executes the compiled MEX function with the specified input arguments.

## MEX Files for Generating a Model Time Activity Curve (TAC)

### Inputs:
- `par`: Model parameters.
- `scant`: Scan time data.
- `blood`: Blood input data.
- `wblood`: Whole blood data.
- `dk`: Decay constant.
- `td`: Time step for numerical convolution calculation.

### Outputs:
- `c`: Computed time activity curve (TAC).
- `s`: Jacobian matrix (if requested).

### MEX Files

1. **ktac_1tcm_mex.cpp**
   - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a one-tissue compartmental model (1TCM).

   - **Compilation Instruction**:
     ```matlab
     mex ktac_1tcm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_1tcm
     ```
   - **Usage**:
     ```matlab
     ktac_1tcm(par, scant, blood, wblood, dk, td)
     ```

2. **ktac_2tcm_mex.cpp**
   - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a two-tissue compartmental model (2TCM).

   - **Compilation Instruction**:
     ```matlab
     mex ktac_2tcm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_2tcm
     ```
   - **Usage**:
     ```matlab
     ktac_2tcm(par, scant, blood, wblood, dk, td)
     ```

3. **ktac_liver_mex.cpp**
    - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a liver kinetic model.

    - **Compilation Instruction**:
      ```matlab
      mex ktac_liver_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_liver
      ```
    - **Usage**:
      ```matlab
      ktac_liver(par, scant, blood, wblood, dk, td)
      ```

4. **ktac_srtm_mex.cpp**
    - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a Simplified Reference Tissue Model (SRTM).

    - **Compilation Instruction**:
      ```matlab
      mex ktac_srtm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_srtm
      ```
    - **Usage**:
      ```matlab
      ktac_srtm(par, scant, blood, wblood, dk, td)
      ```

## MEX Files for TAC Fitting

### Inputs:
- `tac`: Time activity curve (TAC) data.
- `w`: Weights for fitting the TAC data.
- `scant`: Scan time data.
- `blood`: Blood input function.
- `wblood`: Whole blood data.
- `dk`: Decay constant.
- `pinit`: Initial parameters for the model.
- `lb`: Lower bounds for the parameters.
- `ub`: Upper bounds for the parameters.
- `psens`: Sensitivity matrix for the parameters.
- `maxit`: Maximum number of iterations.
- `td`: Time step for numerical convolution calculation.

### Outputs:
- `p`: Estimated parameters.
- `c`: Fitted curve.

### MEX Files

1. **kfit_1tcm_mex.cpp**
   - **Purpose**: Implements the fitting of a one-tissue compartmental model (1TCM) using the Levenberg-Marquardt algorithm.

   - **Compilation Instruction**:
     ```matlab
     mex kfit_1tcm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_1tcm
     ```
   - **Usage**:
     ```matlab
     kfit_1tcm(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

2. **kfit_2tcm_mex.cpp**
   - **Purpose**: Implements the fitting of a two-tissue compartmental model (2TCM) using the Levenberg-Marquardt algorithm for running in MATLAB.

   - **Compilation Instruction**:
     ```matlab
     mex kfit_2tcm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2tcm
     ```
   - **Usage**:
     ```matlab
     kfit_2tcm(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td)
     ```

3. **kfit_liver_mex.cpp**
   - **Purpose**: Implements the fitting of a liver kinetic model using the Levenberg-Marquardt algorithm for running in MATLAB.

   - **Compilation Instruction**:
     ```matlab
     mex kfit_liver_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_liver
     ```
   - **Usage**:
     ```matlab
     kfit_liver(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td)
     ```
     
4. **kfit_srtm_mex.cpp**
   - **Purpose**: Implements the fitting of a Simplified Reference Tissue Model (SRTM) using the Levenberg-Marquardt algorithm in MATLAB.

   - **Compilation Instruction**:
     ```matlab
     mex kfit_srtm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_srtm
     ```
   - **Usage**:
     ```matlab
     kfit_srtm(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td)
     ```

## **MEX Files for Parallel Computing**

## Why OpenMP (OMP) is Used

OpenMP (OMP) is also employed in several MEX files to leverage parallel processing with multiple threads, significantly speeding up the time activity curve fitting process when working with a large number of image voxels. By default, the OMP-enabled MEX files are used for better performance. However, non-OMP versions are also available for users who prefer or require single-threaded execution, as provided above.

## MEX Files for TAC Fitting with OpenMP Parallel Processing

### Inputs:
- `tac`: Time activity curve (TAC) data.
- `w`: Weights for fitting the TAC data.
- `scant`: Scan time data.
- `blood`: Blood input function.
- `wblood`: Whole blood data.
- `dk`: Decay constant.
- `pinit`: Initial parameters for the model.
- `lb`: Lower bounds for the parameters.
- `ub`: Upper bounds for the parameters.
- `psens`: Sensitivity matrix for the parameters.
- `maxit`: Maximum number of iterations.
- `td`: Time step for numerical convolution calculation.

### Outputs:
- `p`: Estimated parameters.
- `c`: Fitted curve.

### MEX Files

1. **kfit_1tcm_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a one-tissue compartmental model (1TCM) using the Levenberg-Marquardt algorithm with OpenMP for running in MATLAB.

   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_1tcm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_1tcm_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_1tcm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_1tcm_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```
   - **Usage**:
     ```matlab
     kfit_1tcm_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

2. **kfit_2tcm_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a two-tissue kinetic model (2TCM) using the Levenberg-Marquardt algorithm with OpenMP for parallel processing in MATLAB.
   
   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_2tcm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2tcm_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_2tcm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2tcm_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```
   - **Usage**:
     ```matlab
     kfit_2tcm_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

3. **kfit_liver_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a liver kinetic model using the Levenberg-Marquardt algorithm with OpenMP parallelization in MATLAB.

   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_liver_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_liver_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_liver_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_liver_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```
   - **Usage**:
     ```matlab
     kfit_liver_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

4. **kfit_srtm_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a Simplified Reference Tissue Model (SRTM) using the Levenberg-Marquardt algorithm with OpenMP for parallel processing in MATLAB.

   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_srtm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_srtm_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_srtm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_srtm_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```


