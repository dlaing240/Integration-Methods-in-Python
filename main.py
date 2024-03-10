import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
# Imports from the other files
from integral import Integral
from one_dim_functions import func_list_1d, func_str_list_1d
from n_dim_functions import func_list_nd, func_desc_nd
from time_tests import time_div_dependence, time_test
from accuracy_test import error_div_dependence, error_div_dependence_nd


#############
#   SETUP   #
#############

# Example variables to use for testing integration methods. These can be changed
x_min = 0
x_max = 10

# Subdivision numbers which can be used to test how accuracy or time varies with number of division. (Can be changed)
divs = np.arange(1, 80, 1)
divs_acc = np.arange(1, 10, 1)
divs_nd = np.arange(1, 10, 1)  # If this is large then higher dimension integrals will take a long time to calculate
divs_nd_acc = np.arange(1, 10, 1)
samples = np.arange(1, 200, 1)

# Initialise Integrator class for each of the 1d functions
integrators_1d = []
for func in func_list_1d:
    integrators_1d.append(Integral(func, a=x_min, b=x_max))

# Initialise Integrator class for each of the nd functions.
integrators_nd = []
# elements of func_list_nd are tuples containing the function and the number of dimensions of the function
for func, d in func_list_nd:
    # create the lower and upper bounds
    a = np.zeros(d)
    b = np.full(d, 1)
    integrators_nd.append(Integral(func, a, b))

# func_nums gives the indices of the above integrator instances to use for each test.
# Each test procedure corresponds to an element in func_nums. Empty elements means the corresponding tests are skipped.
func_nums = [
    [0],  # 1d time test (simpsons and midpoint)
    [4],  # nd time test (simpsons and midpoint)
    [],  # 1d time test (Monte Carlo) - empty because it is skipped by default
    [5],  # nd time test (Monte Carlo)
    [1, 2],  # 1d error test (simpsons and midpoint)
    [],  # nd error test (simpsons and midpoint)
    [5],  # nd error test (Monte Carlo)
    [0, 1, 2, 3, 4, 5]  # Adaptive methods test
]

# Use this instead to test only the adaptive methods (no plots)
# func_nums = [[], [], [], [], [], [], [], [0, 1, 2, 3, 4, 5]]

##################
#   TIME TESTS   #
##################

test_figures = []

# 1D composite midpoint and Simpson's rule:
figures_1d_time_test = []
for i in func_nums[0]:
    cm_durations = time_div_dependence(integrator_method=integrators_1d[i].composite_midpoint, division_nums=divs, repeats=100,
                                       plot=False)
    simpsons_durations = time_div_dependence(integrator_method=integrators_1d[i].simpsons, division_nums=divs,
                                             repeats=100, plot=False)

    fig_1d = plt.figure()
    plt.plot(divs, cm_durations, label="Composite Midpoint")
    plt.plot(divs, simpsons_durations, label="Simpson's Rule")
    plt.legend()
    plt.xlabel('Number of Subdivisions')
    plt.ylabel('Completion Time (s)')
    plt.title(f'Time Performance Comparison, Integrating {func_str_list_1d[i]}')
    figures_1d_time_test.append(fig_1d)

# Time performance of the n_dim integrators:
figures_nd_time_test = []
for i in func_nums[1]:
    cm_nd_durations = time_div_dependence(integrators_nd[i].composite_midpoint_n_dim, division_nums=divs_nd, repeats=1, plot=False)
    simpsons_nd_durations = time_div_dependence(integrators_nd[i].simpsons_n_dim, division_nums=divs_nd, repeats=1,  plot=False)

    fig_nd = plt.figure()
    plt.plot(divs_nd, cm_nd_durations, label="Composite Midpoint")
    plt.plot(divs_nd, simpsons_nd_durations, label="Simpson's")
    plt.title(f'Time Performance Comparison, Integrating {func_desc_nd[i]}')
    plt.legend()
    plt.xlabel('Number of Subdivisions')
    plt.ylabel('Completion Time (s)')
    figures_nd_time_test.append(fig_nd)

# 1d Monte Carlo
mc_fig_list = []
for i in func_nums[2]:
    mc_durations = time_div_dependence(integrator_method=integrators_1d[i].monte_carlo, division_nums=samples, plot=False)
    fig = plt.figure()
    plt.plot(samples, mc_durations, label="Monte Carlo")
    plt.legend()
    plt.xlabel('Number of samples')
    plt.ylabel('Completion time (s)')
    plt.title(f'Completion time vs number of samples integrating {func_str_list_1d[i]} using the Monte Carlo method')
    mc_fig_list.append(fig)

# nd Monte Carlo
mc_fig_list_nd = []
for i in func_nums[3]:
    mc_durations = time_div_dependence(integrator_method=integrators_nd[i].monte_carlo_n_dim, division_nums=samples, plot=False)
    mc_errors = error_div_dependence_nd(integrator=integrators_nd[i], integrator_method=integrators_nd[i].monte_carlo_n_dim,
                                        divs=samples, plot=False)
    fig = plt.figure()
    plt.plot(samples, mc_durations, label="Monte Carlo")
    plt.legend()
    plt.xlabel('Number of Samples')
    plt.ylabel('Completion time (s)')
    plt.title(f'Completion Time vs Number of Samples, Integrating {func_desc_nd[i]} Using the Monte Carlo Method')
    mc_fig_list_nd.append(fig)

###################
#   Error tests   #
###################

# 1d midpoint and simpson's
figs_acc = []
for i in func_nums[4]:
    cm_error = error_div_dependence(integrators_1d[i], integrators_1d[i].composite_midpoint, divs_acc, plot=False)
    sim_error = error_div_dependence(integrators_1d[i], integrators_1d[i].simpsons, divs_acc, plot=False)
    fig = plt.figure()
    plt.plot(divs_acc, cm_error, label='Composite Midpoint')
    plt.plot(divs_acc, sim_error, label='Simpsons')
    plt.legend()
    plt.xlabel('Number of Subdivisions')
    plt.ylabel('Fractional Error')
    plt.title(f'Error vs Number of Subdivisions, Integrating {func_str_list_1d[i]}')
    figs_acc.append(fig)

# nd midpoint and simpson's
for i in func_nums[5]:
    cm_nd_error = error_div_dependence_nd(integrators_nd[i], integrators_nd[i].composite_midpoint_n_dim, divs_nd_acc, plot=False)
    sim_nd_error = error_div_dependence_nd(integrators_nd[i], integrators_nd[i].simpsons_n_dim, divs_nd_acc, plot=False)
    fig = plt.figure()
    plt.plot(divs_nd_acc, cm_nd_error, label='Composite Midpoint')
    plt.plot(divs_nd_acc, sim_nd_error, label='Simpsons')
    plt.legend()
    plt.xlabel('Number of Subdivisions')
    plt.ylabel('Fractional Error')
    plt.title(f'Error vs subdivision number comparison, integrating {func_desc_nd[i]}')
    figs_acc.append(fig)

# nd Monte Carlo
mc_fig_list_error_nd = []
for i in func_nums[6]:
    mc_errors = error_div_dependence_nd(integrator=integrators_nd[i], integrator_method=integrators_nd[i].monte_carlo_n_dim,
                                        divs=samples, plot=False)
    fig = plt.figure()
    plt.plot(samples, mc_errors, label="Monte Carlo")
    plt.legend()
    plt.xlabel('Number of Samples')
    plt.ylabel('Fractional error')
    plt.title(f'Error vs Number of Samples Integrating {func_desc_nd[i]} Using the Monte Carlo Method')
    mc_fig_list_error_nd.append(fig)

#############################
#   Test adaptive methods   #
#############################

print("Testing adaptive methods")
for i in func_nums[7]:
    # Check result using scipy.integrate.quad. Function and lower and upper bounds are given by the instance attributes.
    simpsons_test = time_test(integrators_1d[i].simpsons_adaptive, 0.01)
    midpoint_test = time_test(integrators_1d[i].composite_midpoint_adaptive, 0.01)
    print(
        "\nfunction: ", func_str_list_1d[i],
        "\nExpected Result: ", integrate.quad(integrators_1d[i].func, integrators_1d[i].a[0], integrators_1d[i].b[0])[0],
        "\nAdaptive Simpsons using an error threshold of 1%:",
        "\nResult: ", simpsons_test[1][0], "Completion time: ", simpsons_test[0], "Subdivisions: ", simpsons_test[1][1],
        "\nAdaptive Composite Midpoint using an error threshold of 1%:",
        "\nResult: ", midpoint_test[1][0], "Completion time: ", midpoint_test[0], "Subdivisions: ", midpoint_test[1][1],
          )

# Plot all the results
plt.show()
