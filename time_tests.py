import matplotlib.pyplot as plt
import time


# Functions to demonstrate time performance
def time_test(integrator_method, *method_params):
    """
    Function to measure completion time of an integrator method.

    Parameters
    ----------
    integrator_method:
        Use instance.method for an instance of the integrator class and the desired method.
    method_params:
        Parameters required by the method (if any).

    Returns
    -------
    A tuple of the completion time in seconds and the integral result.

    """
    start = time.perf_counter()
    result = integrator_method(*method_params)
    stop = time.perf_counter()
    return stop - start, result


# integrator_class, func, x_min, x_max,
def time_div_dependence(integrator_method, division_nums, repeats=1, plot=False):
    """
    Function to test the dependence of completion time on the number of subdivisions for an integration method.
    Can be used with non-adaptive methods.
    In the case of Monte Carlo, division_nums instead refers to the number of random samples.

    Parameters
    ----------

    integrator_method:
        Use instance.method for an instance of the integrator class and the desired method.
    division_nums:
        Array of the subdivision numbers to measure completion time with.
    repeats
        How many times the time test repeats for each subdivision.
    plot:
        True/False. Choose whether to produce a plot of the results.

    Returns
    -------
    Array of the completion times.
    """
    durations = []
    for n in division_nums:
        temp_durations = []
        for _ in range(repeats):
            temp_durations.append(time_test(integrator_method, n)[0])
        mean_duration = sum(temp_durations) / repeats
        durations.append(mean_duration)
        # durations.append(time_test(integrator_method, n)[0])
    if plot:
        plt.plot(division_nums, durations)
        plt.show()
    return durations
