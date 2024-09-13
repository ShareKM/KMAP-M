# Kinetic Modeling and Analysis Package (KMAP)

## Overview of the KMAP Library

The KMAP (Kinetic Modeling and Analysis Package) library provides a suite of code for implementing and applying various kinetic models to analyze dynamic positron emission tomography (PET) data. This library includes C/C++ code and MEX files that facilitate seamless integration with MATLAB, enabling efficient kinetic analysis in a flexible and user-friendly environment.

### Supported Kinetic Models

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

### Optimization Methods Used

The KMAP library employs several advanced optimization techniques to accurately fit kinetic models to time-activity curve (TAC) data derived from PET imaging.

1. **Levenberg-Marquardt Algorithm**:
   - **Purpose**: The Levenberg-Marquardt (LM) algorithm is a widely used classic optimization method for solving nonlinear least squares problems, particularly in the context of time activity curve (TAC) fitting. In kinetic modeling, it plays a crucial role in estimating the model parameters so that the computed TAC closely matches the measured TAC data. The algorithm effectively bridges the gap between the steepest descent method and the Gauss-Newton algorithm, providing a balanced approach that is both robust and efficient.
   
   - **How It Works**: 
     - The LM algorithm iteratively updates the parameters by minimizing the sum of the squared differences between the observed and model-predicted TAC values. It combines the advantages of two methods: 
       - **Gradient Descent**: Useful when the parameters are far from their optimal values, where it uses the gradient of the cost function to guide the search direction.
       - **Gauss-Newton**: More effective when the parameters are close to the optimal values, using a second-order approximation of the cost function to make more precise adjustments.
     - The LM algorithm introduces a damping factor that adjusts the step size based on the current position in the parameter space. When far from the minimum, the algorithm behaves more like a gradient descent (larger damping), and as it approaches the minimum, it shifts towards the Gauss-Newton method (smaller damping).

   - **Implementation**: 
     - In `kinlib_optimization.cpp`, the LM algorithm is implemented through the `kmap_levmar` function. This function takes in the initial parameter estimates and iteratively refines them by minimizing the residuals (differences between measured and predicted TACs). 
     - To handle the complexity of multiple parameter constraints, the algorithm also uses a bounded coordinate descent method for solving the intermediate quadratic optimization problem, which ensures that parameter updates stay within specified bounds, improving the stability and reliability of the optimization process.
     - The `kmap_levmar` function also incorporates mechanisms to dynamically adjust the damping factor, improving convergence speed and accuracy. The function iteratively recalculates the Jacobian matrix, which represents the sensitivity of the TAC to each model parameter, ensuring that the parameter updates are optimally directed.

2. **Bounded Coordinate Descent**:
   - **Purpose**: This method uses a coordinate descent optimization strategy to solve the intermediate quadratic optimization problem during the LM optimization. It is particularly effective in situations where the parameters are subject to specific constraints (bounds). 
   
   - **Implementation**: In `kinlib_optimization.cpp`, the `BoundQuadCD` function performs this optimization, iterating over each parameter to minimize the quadratic cost function under the constraints of the parameter bounds. This ensures that the optimization process remains stable and that the parameter estimates are physically meaningful.

3. **Convolution of an Exponential Function and Input Function**:
   - **Purpose**: Convolution of a single exponential function and the blood input function is a building block for the calculation of the analytical solution for different compartmental models.
   
   - **Implementation**: The `kconv_exp` function in `kinlib_common.cpp` efficiently handles the numerical calculation of the convolution for computing the model TACs and their sensitivities concerning different model parameters. This operation directly influences the accuracy of the model fitting process, as it underpins the evaluation of how well the model's predictions align with the observed data.

These optimization and modeling techniques are integral to the KMAP library's ability to accurately model time activity curves in dynamic PET. By employing a combination of the Levenberg-Marquardt algorithm, bounded coordinate descent, and convolution operations, the library ensures robust and reliable parameter estimation, facilitating meaningful insights into tracer kinetics.

## Repository Structure

### 1. **Header File: `kinlib.h`**
   - **Purpose**: Contains function declarations essential for the kinetic modeling process.
   - **Main Functions**:
     - `tac_eval`: Declares the function for evaluating time-activity curves (TACs).
     - `jac_eval`: Declares the function for calculating Jacobian matrices (sensitivity of TACs to model parameters).
     - `frame`: Declares the function for averaging activity over specified time frames.
     - `kconv_exp`: Declares the function for exponential convolution operations.
     - `kmap_levmar`: Declares the Levenberg-Marquardt algorithm for non-linear least squares optimization.

### 2. **Source Files**:
   - **`kinlib_models.cpp`**: 
     - **Purpose**: Contains the implementations of TACs and Jacobian calculations for the various kinetic models.
     - **Key Functions**:
       - **TAC Evaluation**:
         - `tac_eval`: Implements the computation of TACs for the given model parameters.
         - `kconv_1tcm_tac`, `kconv_2tcm_tac`, `kconv_srtm_tac`, `kconv_liver_tac`: Specialized TAC computation functions for the respective models.
       - **Jacobian Calculation**:
         - `jac_eval`: Implements the computation of the Jacobian matrix.
         - `kconv_1tcm_jac`, `kconv_2tcm_jac`, `kconv_srtm_jac`, `kconv_liver_jac`: Specialized Jacobian calculation functions for the respective models.

   - **`kinlib_optimization.cpp`**:
     - **Purpose**: Contains the implementations of optimization algorithms.
     - **Key Functions**:
       - **Optimization**:
         - `kmap_levmar`: Implements the Levenberg-Marquardt algorithm to fit a kinetic model to the TAC data by minimizing the difference between measured and modeled TACs.
         - `BoundQuadCD`: Implements the bounded coordinate descent method for quadratic optimization under constraints.

   - **`kinlib_common.cpp`**:
     - **Purpose**: Contains common functions used across different models and optimizations.
     - **Key Functions**:
       - **Convolution Operations**:
         - `kconv_exp`: Performs the convolution of the input function with a single exponential function, a key building block for more complex compartmental models.
       - **Time Delay Correction**:
         - `time_delay_tac`: Computes time-delayed TAC curves by shifting the TAC in time and interpolating values based on the specified delay.
         - `time_delay_jac`: Computes the Jacobian (gradient) for time delay correction, allowing for the calculation of sensitivity with respect to time delay.
       - **Other Utility Functions**:
         - `frame`: Computes the average activity within specified time frames.
         - `vecnorm2`, `vecnormw`: Helper functions for vector norm calculations.
       - **Parameter Handling**:
         - `setkin`, `getkin`: Functions for defining and handling sensitive parameters during optimization.  

       

### 3. **MEX Files for MATLAB**

This directory contains several MEX files for performing kinetic modeling in MATLAB. Precompiled MEX binaries (`*.mexa64` files) are also included and can be used directly in MATLAB without recompilation, provided they are compatible with your system architecture.

Each MEX file is accompanied by detailed instructions on its usage, input, output, and compilation. Please refer to the comments within each MEX source file (`*.cpp`) for precise information on how to work with these files.

### Compilation Instruction:

The compilation command should be run in the MATLAB Command Window or a terminal/shell where you have access to the MATLAB mex compiler. This command compiles the C++ source code into a MEX file that can be executed within MATLAB.

### Usage:

The usage command should be run within the MATLAB environment after the MEX file has been successfully compiled. It executes the compiled MEX function with the specified input arguments.

### MEX Files for Generating a Model Time Activity Curve (TAC)

0. **Input and Outputs**:
   - **Inputs**:
     - `par`: Model parameters.
     - `scant`: Scan time data.
     - `blood`: Blood input data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `td`: Time step for numerical convolution calculation.
   - **Outputs**: 
     - `c`: Computed time activity curve (TAC).
     - `s`: Jacobian matrix (if requested).
       
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

### MEX Files for TAC Fitting
0. **Inputs and Outputs**
   All the input and output variables are the same for different models
   - **Inputs**: 
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
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time step for numerical convolution calculation.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.
       
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

### 4. **MEX Files for Parallel Computing**

### Why OpenMP (OMP) is Used

OpenMP (OMP) is also employed in several MEX files to leverage parallel processing with multiple threads, significantly speeding up the time activity curve fitting process when working with a large number of image voxels. By default, the OMP-enabled MEX files are used for better performance. However, non-OMP versions are also available for users who prefer or require single-threaded execution, as provided above.

### MEX Files for TAC Fitting with OpenMP Parallel Processing
0. **Inputs and Outputs**
   All the input and output variables are the same for different models
   - **Inputs**: 
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
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time step for numerical convolution calculation.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.

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
