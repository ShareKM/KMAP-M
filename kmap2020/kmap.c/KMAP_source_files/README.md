# KMAP Source Files

This folder contains the core C++ source code files used to define and implement the kinetic models, optimization routines, and utility functions. These files are compiled and integrated with the MEX functions to allow for MATLAB integration.

## Source Files

1. **`kinlib_models.cpp`**: 
   - **Purpose**: Contains the implementations of TAC (Time-Activity Curve) and Jacobian calculations for the various kinetic models.
   - **Key Functions**:
     - **TAC Evaluation**:
       - `tac_eval`: Computes TACs for given model parameters.
       - `kconv_1tcm_tac`, `kconv_2tcm_tac`, `kconv_srtm_tac`, `kconv_liver_tac`: Specialized TAC computation functions for the respective models.
     - **Jacobian Calculation**:
       - `jac_eval`: Computes the Jacobian matrix for TACs.
       - `kconv_1tcm_jac`, `kconv_2tcm_jac`, `kconv_srtm_jac`, `kconv_liver_jac`: Specialized Jacobian calculation functions.

2. **`kinlib_optimization.cpp`**:
   - **Purpose**: Contains the optimization routines used in the library, including the Levenberg-Marquardt algorithm.
   - **Key Functions**:
     - **Optimization**:
       - `kmap_levmar`: Levenberg-Marquardt algorithm for optimizing model parameters to fit TAC data.
       - `BoundQuadCD`: Bounded coordinate descent for quadratic optimization under constraints.

3. **`kinlib_common.cpp`**:
   - **Purpose**: Contains common utility functions used across different models and optimization routines.
   - **Key Functions**:
     - **Convolution Operations**:
       - `kconv_exp`: Convolution of an input function with a single exponential, a key component in compartmental models.
     - **Time Delay Correction**:
       - `time_delay_tac`: Computes time-delayed TAC curves by shifting the TAC in time.
       - `time_delay_jac`: Computes the Jacobian for time delay correction.
     - **Other Utility Functions**:
       - `frame`: Computes average activity within specified time frames.
       - `vecnorm2`, `vecnormw`: Helper functions for vector norm calculations.
     - **Parameter Handling**:
       - `setkin`, `getkin`: Functions for defining and handling sensitive parameters during optimization.

## Optimization Methods Used

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

## Compilation

To compile these files for integration with MATLAB, use the MATLAB `mex` command from within the `KMAP_mex_files` folder.

Refer to the detailed compilation instructions in `KMAP_mex_files/README.md`

