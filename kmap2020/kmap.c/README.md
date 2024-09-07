# Kinetic Modeling Library (KMAP)

## Introduction to Kinetic Models

Kinetic models are mathematical frameworks used to describe the dynamic behavior of biological systems, particularly the movement and transformation of substances within the body over time. In medical imaging, such as Positron Emission Tomography (PET), kinetic models are essential for interpreting the distribution and concentration of radiotracers. These tracers are introduced into the body and tracked as they move through different compartments (e.g., blood, tissues) and undergo biochemical reactions.

### Why We Need Kinetic Models

1. **Quantification of Biological Processes**: Kinetic models allow for the quantification of biological processes such as metabolism, blood flow, and receptor-ligand binding, which are crucial for diagnosing and monitoring diseases.
  
2. **Non-invasive Insight**: Through imaging techniques like PET, kinetic models provide non-invasive insights into the functioning of organs and tissues, helping in early detection and accurate monitoring of diseases like cancer, Alzheimer's, and cardiovascular disorders.

3. **Understanding Tracer Kinetics**: By fitting kinetic models to time-activity curves (TACs) derived from PET data, researchers can extract parameters that describe the rates of tracer uptake, binding, and clearance. These parameters are invaluable for clinical and research applications.

## Overview of the KMAP Library

The KMAP (Kinetic Modeling Analysis Package) library provides a suite of tools for implementing and applying various kinetic models to PET data. This library includes C++ code and MEX files that facilitate seamless integration with MATLAB, enabling efficient kinetic analysis in a flexible and user-friendly environment.

### Supported Kinetic Models

1. **1-Tissue 3-Parameter Model (1T3P)**
   - **Description**: A simplified model where the tracer exchanges between plasma and a single tissue compartment.
   - **Parameters**: \( K_1 \) (influx rate), \( k_2 \) (efflux rate), \( V_b \) (vascular fraction).

2. **2-Tissue 5-Parameter Model (2T5P)**
   - **Description**: A more detailed model with two compartments representing different tissue types or states, providing a richer description of tracer kinetics.
   - **Parameters**: \( K_1 \), \( k_2 \), \( k_3 \) (binding rate), \( k_4 \) (dissociation rate), \( V_b \).

3. **Simplified Reference Tissue Model (SRTM)**
   - **Description**: This model is used when a reference region devoid of the receptor of interest is available, helping to estimate binding potential without requiring blood samples.
   - **Parameters**: \( R_1 \) (relative delivery rate), \( k_2 \), \( BP_{nd} \) (binding potential), \( V_b \).

4. **Liver Model**
   - **Description**: A specialized model accounting for the dynamics in liver tissue, where both plasma and whole blood compartments are considered.
   - **Parameters**: \( K_1 \), \( k_2 \), \( k_3 \), \( k_4 \), \( K_a \), \( f_a \) (fractional arterial blood), \( V_b \).

### Optimization Methods Used

The KMAP library employs several advanced optimization techniques to accurately fit kinetic models to time-activity curve (TAC) data derived from PET imaging.

1. **Levenberg-Marquardt Algorithm**:
   - **Purpose**: The Levenberg-Marquardt (LM) algorithm is a widely used method for solving nonlinear least squares problems, particularly in the context of curve fitting. In kinetic modeling, it plays a crucial role in adjusting the model parameters so that the computed TAC closely matches the observed data. The algorithm effectively bridges the gap between the steepest descent method and the Gauss-Newton algorithm, providing a balanced approach that is both robust and efficient.
   
   - **How It Works**: 
     - The LM algorithm iteratively updates the parameters by minimizing the sum of the squared differences between the observed and model-predicted TAC values. It combines the advantages of two methods: 
       - **Gradient Descent**: Useful when the parameters are far from their optimal values, where it uses the gradient of the cost function to guide the search direction.
       - **Gauss-Newton**: More effective when the parameters are close to the optimal values, using a second-order approximation of the cost function to make more precise adjustments.
     - The LM algorithm introduces a damping factor that adjusts the step size based on the current position in the parameter space. When far from the minimum, the algorithm behaves more like a gradient descent (larger damping), and as it approaches the minimum, it shifts towards the Gauss-Newton method (smaller damping).

   - **Implementation**: 
     - In `kinlib_optimization.cpp`, the LM algorithm is implemented through the `kmap_levmar` function. This function takes in the initial parameter estimates and iteratively refines them by minimizing the residuals (differences between observed and predicted TACs). 
     - To handle the complexity of multiple parameters, the algorithm also uses a bounded coordinate descent method, which ensures that parameter updates stay within specified bounds, improving the stability and reliability of the optimization process.
     - The `kmap_levmar` function also incorporates mechanisms to dynamically adjust the damping factor, improving convergence speed and accuracy. The function iteratively recalculates the Jacobian matrix, which represents the sensitivity of the TAC to each model parameter, ensuring that the parameter updates are optimally directed.

2. **Bounded Coordinate Descent**:
   - **Purpose**: This method is employed to optimize individual parameters while keeping the others fixed. It is particularly effective in situations where the parameters are subject to specific constraints (bounds). The method complements the LM algorithm by refining the parameter estimates within the allowable range.
   
   - **Implementation**: In `kinlib_optimization.cpp`, the `BoundQuadCD` function performs this optimization, iterating over each parameter to minimize the cost function under the constraints of the parameter bounds. This ensures that the optimization process remains stable and that the parameter estimates are physically meaningful.

3. **Convolution Operations**:
   - **Purpose**: Convolution with exponential functions is essential in kinetic modeling to simulate the tracer's behavior as it transitions between different biological compartments over time.
   
   - **Implementation**: The `kconv_exp` function in `kinlib_common.cpp` efficiently handles these convolutions, which are integral to computing the TACs and their sensitivities. This operation directly influences the accuracy of the model fitting process, as it underpins the evaluation of how well the model's predictions align with the observed data.

These optimization techniques are integral to the KMAP library's ability to accurately model biological processes observed in PET imaging. By employing a combination of the Levenberg-Marquardt algorithm, bounded coordinate descent, and convolution operations, the library ensures robust and reliable parameter estimation, facilitating meaningful insights into tracer kinetics.

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
     - **Purpose**: Contains the implementations of TAC and Jacobian calculations for the various kinetic models.
     - **Key Functions**:
       - **TAC Evaluation**:
         - `tac_eval`: Implements the computation of TACs for the given model parameters.
         - `kconv_1t3p_tac`, `kconv_2t5p_tac`, `kconv_srtm_tac`, `kconv_liver_tac`: Specialized TAC computation functions for the respective models.
       - **Jacobian Calculation**:
         - `jac_eval`: Implements the computation of the Jacobian matrix.
         - `kconv_1t3p_jac`, `kconv_2t5p_jac`, `kconv_srtm_jac`, `kconv_liver_jac`: Specialized Jacobian calculation functions for the respective models.

   - **`kinlib_optimization.cpp`**:
     - **Purpose**: Contains the implementations of optimization algorithms.
     - **Key Functions**:
       - **Optimization**:
         - `kmap_levmar`: Implements the Levenberg-Marquardt algorithm to fit the kinetic model to the TAC data by minimizing the difference between observed and modeled TACs.
         - `BoundQuadCD`: Implements the bounded coordinate descent method for parameter optimization under constraints.
         - `lema_gsn`: Implements the Levenberg-Marquardt algorithm with additional constraints.

   - **`kinlib_common.cpp`**:
     - **Purpose**: Contains common functions used across different models and optimizations.
     - **Key Functions**:
       - **Convolution Operations**:
         - `kconv_exp`: Performs the convolution of input functions with an exponential decay, a key operation in kinetic modeling.
       - **Other Utility Functions**:
         - `frame`: Computes the average activity within specified time frames.
         - `vecnorm2`, `vecnormw`: Helper functions for vector norm calculations.
       - **Parameter Handling**:
         - `setkin`, `getkin`: Functions for handling sensitive parameters during optimization.  
       

### 3. **MEX Files for each of the models**

## MEX Files

This directory contains several MEX files for performing kinetic modeling in MATLAB. Precompiled MEX binaries (`*.mexa64` files) are also included and can be used directly in MATLAB without recompilation, provided they are compatible with your system architecture.

Each MEX file is accompanied by detailed instructions on its usage, input, output, and compilation. Please refer to the comments within each MEX source file (`*.cpp`) for precise information on how to work with these files.

### Why OpenMP (OMP) is Used

OpenMP (OMP) is employed in several MEX files to leverage parallel processing, significantly speeding up the fitting process when working with large datasets or multiple voxels. By default, the OMP-enabled MEX files are used for better performance. However, non-OMP versions are also available for users who prefer or require single-threaded execution.

### Compilation Instruction:

The compilation command should be run in the MATLAB Command Window or a terminal/shell where you have access to the MATLAB mex compiler. This command compiles the C++ source code into a MEX file that can be executed within MATLAB.

### Usage:

The usage command should be run within the MATLAB environment after the MEX file has been successfully compiled. It executes the compiled MEX function with the specified input arguments.

### MEX Files Overview

1. **kfit_1t3p_mex.cpp**
   - **Purpose**: Implements the fitting of a one-tissue kinetic model (1T3P) using the Levenberg-Marquardt algorithm in MATLAB.
   - **Inputs**: 
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `lb`: Lower bounds for the parameters.
     - `ub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.
   - **Compilation Instruction**:
     ```matlab
     mex kfit_1t3p_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_1t3p
     ```
   - **Usage**:
     ```matlab
     kfit_1t3p(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

2. **kfit_1t3p_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a one-tissue kinetic model (1T3P) using the Levenberg-Marquardt algorithm with OpenMP for parallel processing in MATLAB.
   - **Inputs**:
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `lb`: Lower bounds for the parameters.
     - `ub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.

   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_1t3p_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_1t3p_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_1t3p_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_1t3p_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```
   - **Usage**:
     ```matlab
     kfit_1t3p_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

3. **kfit_2t5p_mex.cpp**
   - **Purpose**: Implements the fitting of a two-tissue kinetic model (2T5P) using the Levenberg-Marquardt algorithm in MATLAB.
   - **Inputs**:
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `plb`: Lower bounds for the parameters.
     - `pub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.
   - **Compilation Instruction**:
     ```matlab
     mex kfit_2t5p_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2t5p
     ```
   - **Usage**:
     ```matlab
     kfit_2t5p(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td)
     ```

4. **kfit_2t5p_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a two-tissue kinetic model (2T5P) using the Levenberg-Marquardt algorithm with OpenMP for parallel processing in MATLAB.
   - **Inputs**:
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `lb`: Lower bounds for the parameters.
     - `ub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.
   
   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_2t5p_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2t5p_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_2t5p_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_2t5p_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```
   - **Usage**:
     ```matlab
     kfit_2t5p_mex_omp(tac, w, scant, blood, wblood, dk, pinit, lb, ub, psens, maxit, td)
     ```

5. **kfit_liver_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a liver kinetic model using the Levenberg-Marquardt algorithm with OpenMP parallelization in MATLAB.
   - **Inputs**:
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `lb`: Lower bounds for the parameters.
     - `ub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.

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

6. **kfit_srtm_mex.cpp**
   - **Purpose**: Implements the fitting of a Simplified Reference Tissue Model (SRTM) using the Levenberg-Marquardt algorithm in MATLAB.
   - **Inputs**:
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `plb`: Lower bounds for the parameters.
     - `pub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `p`: Estimated parameters.
     - `c`: Fitted curve.
   - **Compilation Instruction**:
     ```matlab
     mex kfit_srtm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_srtm
     ```
   - **Usage**:
     ```matlab
     kfit_srtm(tac, w, scant, blood, wblood, dk, pinit, plb, pub, psens, maxit, td)
     ```

7. **kfit_srtm_mex_omp.cpp**
   - **Purpose**: Implements the fitting of a Simplified Reference Tissue Model (SRTM) using the Levenberg-Marquardt algorithm with OpenMP for parallel processing in MATLAB.
   - **Inputs**:
     - `tac`: Time activity curve (TAC) data.
     - `w`: Weights for the TAC data.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `pinit`: Initial parameters for the model.
     - `lb`: Lower bounds for the parameters.
     - `ub`: Upper bounds for the parameters.
     - `psens`: Sensitivity matrix for the parameters.
     - `maxit`: Maximum number of iterations for the fitting algorithm.
     - `td`: Time duration for the scan.

   - **Compilation Instructions**:
     - **Windows**:
       ```matlab
       mex kfit_srtm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_srtm_mex_omp CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
       ```
     - **Linux**:
       ```matlab
       mex kfit_srtm_mex_omp.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output kfit_srtm_mex_omp CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
       ```

8. **ktac_1t3p_mex.cpp**
   - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a one-tissue kinetic model (1T3P) in MATLAB.
   - **Inputs**:
     - `par`: Model parameters.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `c`: Computed time activity curve (TAC).
     - `s`: Jacobian matrix (if requested).
   - **Compilation Instruction**:
     ```matlab
     mex ktac_1t3p_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_1t3p
     ```
   - **Usage**:
     ```matlab
     ktac_1t3p(par, scant, blood, wblood, dk, td)
     ```

9. **ktac_2t5p_mex.cpp**
   - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a two-tissue kinetic model (2T5P) in MATLAB.
   - **Inputs**:
     - `par`: Model parameters.
     - `scant`: Scan time data.
     - `blood`: Blood data.
     - `wblood`: Whole blood data.
     - `dk`: Decay constant.
     - `td`: Time duration for the scan.
   - **Outputs**: 
     - `c`: Computed time activity curve (TAC).
     - `s`: Jacobian matrix (if requested).
   - **Compilation Instruction**:
     ```matlab
     mex ktac_2t5p_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_2t5p
     ```
   - **Usage**:
     ```matlab
     ktac_2t5p(par, scant, blood, wblood, dk, td)
     ```

10. **ktac_liver_mex.cpp**
    - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a liver kinetic model in MATLAB.
    - **Inputs**:
      - `par`: Model parameters.
      - `scant`: Scan time data.
      - `blood`: Blood data.
      - `wblood`: Whole blood data.
      - `dk`: Decay constant.
      - `td`: Time duration for the scan.
    - **Outputs**: 
      - `c`: Computed time activity curve (TAC).
      - `s`: Jacobian matrix (if requested).
    - **Compilation Instruction**:
      ```matlab
      mex ktac_liver_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_liver
      ```
    - **Usage**:
      ```matlab
      ktac_liver(par, scant, blood, wblood, dk, td)
      ```

11. **ktac_srtm_mex.cpp**
    - **Purpose**: Implements the computation of the time activity curve (TAC) and its Jacobian for a Simplified Reference Tissue Model (SRTM) in MATLAB.
    - **Inputs**:
      - `par`: Model parameters.
      - `scant`: Scan time data.
      - `blood`: Blood data.
      - `wblood`: Whole blood data.
      - `dk`: Decay constant.
      - `td`: Time duration for the scan.
    - **Outputs**: 
      - `c`: Computed time activity curve (TAC).
      - `s`: Jacobian matrix (if requested).
    - **Compilation Instruction**:
      ```matlab
      mex ktac_srtm_mex.cpp kinlib_models.cpp kinlib_optimization.cpp kinlib_common.cpp -output ktac_srtm
      ```
    - **Usage**:
      ```matlab
      ktac_srtm(par, scant, blood, wblood, dk, td)
      ```

