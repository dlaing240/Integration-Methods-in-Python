import scipy.integrate as integrate
import matplotlib.pyplot as plt


def error_div_dependence(integrator, integrator_method, divs, plot=False):
    """
    Function to test the dependence of error on the number of subdivisions for 1D integration methods.

    Parameters
    ----------
    integrator:
        The instance of the Integrator class
    integrator_method:
        Integration method to be tested. Use instance.method
    divs:
        Array of the subdivision numbers to be used.
    plot:
        True/False. Choose whether the results are plotted.

    Returns
    -------
    An array of the errors for each subdivision number tested.
    """
    actual_value = integrate.quad(func=integrator.func, a=integrator.a[0], b=integrator.b[0])[0]
    errors = []
    for n in divs:
        result = integrator_method(n)
        errors.append(abs(result - actual_value) / actual_value)

    if plot:
        plt.plot(divs, errors)
        plt.show()
    return errors


def error_div_dependence_nd(integrator, integrator_method, divs, plot=False):
    """
    Function to test the dependence of error on the number of subdivisions for N-dim integration methods.

    Parameters
    ----------
    integrator:
        The instance of the Integrator class
    integrator_method:
        Integration method to be tested. Use instance.method
    divs:
        Array of the subdivision numbers to be used.
    plot:
        True/False. Choose whether the results are plotted.

    Returns
    -------
    An array of the errors for each subdivision number tested.
    """
    ranges = []
    for i in range(len(integrator.a)):
        ranges.append((integrator.a[i], integrator.b[i]))
    actual_value = integrate.nquad(func=integrator.func, ranges=ranges)

    errors = []
    for n in divs:
        result = integrator_method(n)
        errors.append(abs(result - actual_value[0])/actual_value[0])

    if plot:
        plt.plot(divs, errors)
        plt.show()
    return errors
