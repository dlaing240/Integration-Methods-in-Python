# Integration Methods in Python

## Overview

The focus of this project is to investigate various integration methods and their implementations in Python.
The following methods were investigated:
* Composite Midpoint
* Simpson's 1/3 rule
* Monte Carlo Integration

Each method was also generalised to n dimensions, and adaptive methods for composite midpoint and Simpson's rule
were also made.

The methods are defined in a class called 'Integral'. This class is initialised with a function and the bounds of integration
and the different methods can then be called on the instance to estimate the value of the integral.

Functions to measure the time performance and accuracy are defined in the accuracy_test.py and time_tests.py files. These are
used in the analysis procedure, which is written in the main.py file. There are two lists of functions used for the tests,
one for one-dimensional functions and the other for n-dimensional functions, and each contains functions with a range of 
behaviours to test the integration methods.

### Example Usage

This project requires the following libraries: NumPy, SciPy and Matplotlib.

The integration methods can be used in this way:
```
from integral import Integral
import numpy as np

# Define the function to be integrated
def f(x):
    return np.sin(x)

# Initialise the integral with the function and the bounds of integration
integral = Integral(f, 0, np.pi)

# Use the desired method to estimate the value of the integral
result = integral.composite_midpoint(n=100)  # This uses the composite midpoint method with 100 subdivisions.
```
## Discussion

### Introduction to the Integration Methods

The composite midpoint method is a method of approximating the result of integration by creating a number of rectangles with a fixed width
and calculating their area. This approximation is very simple, but becomes increasingly accurate for higher numbers of rectangles.
The integral calculation for each rectangle involves multiplying the rectangle width by the value of the function evaluated at the midpoint
of that rectangle. These calculations are summed to get the final result of integration.

The Simpson's 1/3 rule is an example of a Newton-Cotes formula, a group of formulas used for numerical integration. It is a more complex calculation
than the midpoint rule, and requires the function to be evaluated at three points instead of just at the midpoint. The approximation is known to be exact for
polynomials of third degree and below. The Simpson's rule can be used by splitting the integration bounds into a number of subdivisions, in the same
way as the midpoint rule, which gives the composite Simpson's method.

Monte Carlo integration is significantly different to the previous two, as it involves randomly sampling values from the uniform distribution.
The function is evaluated at each of these values, and then the mean of the evaluations is multiplied by the full range of the integration to
get the estimate for the integral. The estimate converges towards the true integral as the number of samples increases.

### Generalisation to n dimensions

The midpoint and Simpson's rule can be generalised to n dimensions by establishing an n-dimensional grid, and then evaluating either the midpoint rule
or Simpson's rule for each element of the grid and sum to get the integral estimate. The generalisation of Monte Carlo integration to n dimensions is 
even simpler, as it just involves randomly sampling values for each dimension, and evaluating the function for the resulting n-d coordinate.

### Adaptive Methods

The motivation for adaptive methods is to be able to estimate an integral to a certain desired accuracy efficiently. This requires the optimum
number of subdivisions to be used. For example, if a 2nd order polynomial was being integrated using a composite Simpson's method, a subdivision number
higher than one would result in unnecessary calculations and function evaluations being performed, since Simpson's would be exact for this polynomial.
Adaptive methods therefore recursively use the integration methods until the desired (estimated) accuracy is reached at which point the procedure ends.

The implementation for the adaptive composite midpoint involves recursively calling the composite midpoint method at twice the number of subdivisions.
The composite midpoint only invovles one function evaluation for the midpoint of each subdivision, so this simple implementation
is enough. Simpson's is more expensive however, so it is important to reuse as many calculations and evaulations as possible, so a slightly different implementation
was used.

The accuracy is estimated by checking the difference between the latest calculation and the previous calculation. There are cases where using absolute error
for this purpose is less sensible (such as a high order polynomial), and other cases where using relative error is less sensible (for example, if the difference is
small compared to the value of the integral, but not compared to the actual answer). Relative error is used in these implementations, however an improvement
would be to incorporate both absolute and relative error.

## Analysis

In all cases with the 1D functions, the completion time of integral calculations scales linearly with the numbers of subdivisions, which is unsurprising.
As expected, the composite midpoint method with it's fewer function evaluations has a shorter completion time per subdivision than the composite Simpson's
method. The completion time of Monte Carlo integration depends on the number of samples in a very similar way to the time performance of composite midpoint, as it
uses one function evaluation per sample. However, Monte Carlo is less precise at comparable subdivisions/sample numbers. The exactness of Simpson's is demontsrated
up to polynomials of third degree, as expected, shown by the constant, near zero error (only affected by numerical precision), and this exactness is broken for the
4th degree polynomial.
![simps and midpoint - time - 1d - 20000](https://github.com/dlaing240/Integration-Methods-in-Python/assets/159714200/1f68ac8a-b1dd-4ff3-9c86-e5e8d7449b32)
![simps and midpoint - errors - x3](https://github.com/dlaing240/Integration-Methods-in-Python/assets/159714200/aa6a65ca-8cf7-4275-ad13-f4fb088fa6b9)
![simps and midpoint - errors - x4](https://github.com/dlaing240/Integration-Methods-in-Python/assets/159714200/9072087a-0791-4016-b027-ecea3ffba4e6)



In the n-dimensional case, the completion time for Simpson's and composite midpoint no longer scales linearly, as now calculations must be done in each dimension.
As a result the time performance for these methods suffers for higher dimensions. The n dimensional generalisation of Monte Carlo integration still only uses one function
evaluation per point though, and so even at higher dimensions the completion time is shown to increase linearly with number of samples. This highlights that Monte Carlo
integration has an advantage for estimating integrals that are otherwise especially challenging.
![simps and midpoint - time - 4d](https://github.com/dlaing240/Integration-Methods-in-Python/assets/159714200/5ba701f5-d09d-4583-be00-d12254b17f58)
![mc - time - 5d](https://github.com/dlaing240/Integration-Methods-in-Python/assets/159714200/4a378a2e-a267-41a0-9e4b-c932c04c658c)
![mc - errors - 5d](https://github.com/dlaing240/Integration-Methods-in-Python/assets/159714200/3e09a1b3-6c65-47c5-9f4a-98c15ae0eb82)


The advantage of the adaptive integration methods was also demonstrated in the analysis procedure, using all the 1D test functions with a target error of 1%.
The methods are able to successfully estimate the value of integration to the desired precision, without requiring any prior calculations or intuition as to what
the optimum number of subdivisions would be.

The project demonstrates the effectiveness of each integration method. Simpson's has good accuracy and is even exact for polynomials of third degree and below, while
the composite midpoint method is a simple yet effective method for quick approximations. Monte Carlo integration excels at estimations for functions of higher dimensions,
which are otherwise challenging.


