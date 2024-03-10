import numpy as np
import random


class Integral:
    """
    A class for performing various integration methods.

    Parameters
    ----------
    func:
        The function to be integrated
    a:
        The lower bound of the integral
    b:
        The upper bound of the integral

    Methods
    -------
    composite_midpoint(n):
        Approximates integral using the composite midpoint rule.
    composite_midpoint_adaptive(e, limit):
        Uses the composite midpoint rule adaptively to reach the desired precision.
    composite_midpoint_n_dim(n):
        Uses composite midpoint rule generalised to n dimensions.
    simpsons(n):
        Uses the composite Simpson's rule to estimate the integral.
    simpsons_adaptive(e):
        Adaptively uses Simpsons to estimate the integral with a desired precision efficiently.
    simpsons_n_dim(n):
        Composite Simpsons rule generalised to n dimensions.
    monte_carlo(n):
        Estimates integral using the Monte Carlo method.
    monte_carlo_n_dim(n):
        Monte Carlo method generalised to n dimensions.
    """
    def __init__(self, func, a, b):
        """
        Initialises the integral class for a given function and integration bounds.
        Parameters
        ----------
        func : callable
            The function to be integrated
        a : float or list
            The lower bound(s) for the integration methods. In the case of n-d functions it is a list of the lower
             bounds for each dimension.
        b : float or list
            The upper bound(s) for the integration methods. In the case of n-d functions it is a list of the upper
             bounds for each dimension.
        """
        if not isinstance(a, (list, tuple, np.ndarray)):
            a = [a]
        if not isinstance(b, (list, tuple, np.ndarray)):
            b = [b]

        self.func = func
        self.a = a
        self.b = b
        self.width = None

        if len(a) == 1:
            self.width = b[0] - a[0]  # Used to improve performance of adaptive composite midpoint

        self.simps_sub_div = None
        self.mid_sub_div = None

    def composite_midpoint(self, n):
        """
        Approximates integral using the composite midpoint rule.
        Parameters
        ----------
        n : int
            Number of subdivisions

        Returns
        -------
        Estimate for the integral.
        """

        # division_width = (self.b - self.a) / n  # width for each rectangle
        division_width = self.width / n  # reducing the calculations slightly for adaptive method
        integrals = []
        left = self.a[0]  # the lower bound of the first rectangle
        for i in range(0, n):
            # upper bound of the rectangle is found by adding the division width to the lower bound
            right = left + division_width
            integral = division_width * self.func((left + right) / 2)
            integrals.append(integral)
            left = right  # sets the upper bound to be the lower bound of the next rectangle
        return sum(integrals)

    def composite_midpoint_adaptive(self, e, limit=5000):
        """
        Uses the composite midpoint rule adaptively to reach the desired precision.
        Parameters
        ----------
        e : float
            The relative error tolerance, indicating when to end the integration
        limit : int
            (Optional) The maximum number of subdivisions before the procedure terminates.
            Can use if there's a chance the error threshold isn't reached

        Returns
        -------
        The estimate for the integral.
        """
        def recursive_integration(prev, n):
            curr = self.composite_midpoint(n)
            if abs(prev - curr) / curr < e or n > limit:
                return curr, n
            return recursive_integration(curr, n*2)
        result = recursive_integration(self.composite_midpoint(1), 2)
        return result

    def _establish_grid(self, n):
        """
        A method used by n-dimensional integration methods to set up the n-dimensional grid.
        """
        dimension_number = len(self.a)  # finds number of dimensions of function
        div_widths = []
        start_points = []
        end_points = []
        # find division width for each dimension
        for i in range(dimension_number):
            div_width = (self.b[i] - self.a[i]) / n
            div_widths.append(div_width)
            start_points.append(self.a[i] + div_width / 2)
            end_points.append(self.b[i] - div_width / 2)
        volume = np.prod(div_widths)  # 'volume' of each division element is used in the integral calculation
        return div_widths, volume, start_points, end_points, dimension_number

    def composite_midpoint_n_dim(self, n):
        """
        Uses composite midpoint rule generalised to n dimensions.
        Parameters
        ----------
        n : int
            The number of subdivisions.

        Returns
        -------
        Estimate for the integral
        """
        div_widths, volume, start_points, end_points, dimension_number = self._establish_grid(n)
        current_point = start_points.copy()
        d = 0  # label for the dimension
        integrals = []
        while True:
            # The div_widths[d]/2 term is to mitigate rounding errors.
            while current_point[d] + div_widths[d] >= end_points[d] + div_widths[d]/2:
                d += 1
                if d == dimension_number:
                    integrals.append(volume * self.func(*current_point))
                    total_integral = sum(integrals)
                    return total_integral
            integrals.append(volume * self.func(*current_point))
            current_point[d] += div_widths[d]
            for i in range(0, d):
                # This is to reset the current point to the start of next 'row'.
                current_point[i] = start_points[i]
            d = 0

    def simpsons(self, n=1):
        """
        Uses the composite Simpson's rule to estimate the integral.
        Parameters
        ----------
        n : int
            Number of subdivisions. (Default is 1, if higher, then this is composite simpson's.)

        Returns
        -------
        Estimate for the integral
        """
        division_width = (self.b[0] - self.a[0])/n
        coefficient = division_width / 6
        total_integral = 0
        left = self.a[0]
        term_1 = self.func(left)
        for i in range(0, n):
            right = left + division_width
            term_3 = self.func(right)
            term_2 = 4 * self.func((left + right) / 2)
            integral = coefficient * (term_1 + term_2 + term_3)
            total_integral += integral
            left = right
            term_1 = term_3
        result = total_integral
        return result

    def _simpsons(self, a, fa, b, fb):
        """
        Used by the adaptive Simpson's method to evaluate the Simpson's rule
        """
        m = (a + b)/2
        fm = self.func(m)
        res = abs(b - a) / 6 * (fa + 4 * fm + fb)
        return m, fm, res

    def _simpsons_adaptive_eff(self, a, fa, b, fb, prev, e, m, fm):
        """
        Recursive function for the adaptive Simpson's method.
        """
        # evaluate simpsons for the left subdivision, saving lm & flm for reuse
        lm, flm, left = self._simpsons(a, fa, m, fm)
        # evaluate simpsons for right subdivision, saving rm & frm for reuse
        rm, frm, right = self._simpsons(m, fm, b, fb)
        self.simps_sub_div += 1  # One division becomes two subdivisions.

        if abs(left + right - prev) / (left + right) < e:
            return left + right
        return (self._simpsons_adaptive_eff(a, fa, m, fm, left, e, lm, flm) +
                self._simpsons_adaptive_eff(m, fm, b, fb, right, e, rm, frm))

    def simpsons_adaptive(self, e):
        """
        Adaptively uses Simpsons to estimate the integral with a desired precision efficiently.
        Parameters
        ----------
        e : float
            Relative error tolerance threshold to indicate when the integration procedure terminates.
        Returns
        -------
        Estimate for the integral.
        """
        a, b = self.a[0], self.b[0]
        fa, fb = self.func(a), self.func(b)
        m, fm, res = self._simpsons(a=a, fa=fa, fb=fb, b=b)
        self.simps_sub_div = 1
        return self._simpsons_adaptive_eff(a, fa, b, fb, res, e, m, fm), self.simps_sub_div

    def simpsons_adaptive_ineff(self, cut_off=10**-6, limit=5000):
        """
        The former adaptive simpsons method, found to be less efficient than the replacement.
        """
        def recursive_integration(prev, n):
            curr = self.simpsons(n)
            if np.abs(curr - prev) < cut_off or n > limit:
                return curr, n
            return recursive_integration(curr, n * 2)
        result = recursive_integration(self.simpsons(4), 8)
        return result

    def simpsons_n_dim(self, n):
        """
        Composite Simpsons rule generalised to n dimensions.
        Parameters
        ----------
        n : int
            Number of subdivisions
        Returns
        -------
        Estimate for the integral
        """
        div_widths, volume, start_points, end_points, dimension_number = self._establish_grid(n)
        coefficient = volume/6
        a_start_points = []
        b_start_points = []
        for i in range(dimension_number):
            a_start_points.append(self.a[i])
            b_start_points.append(self.a[i] + div_widths[i])
        current_point = start_points.copy()
        current_a = a_start_points.copy()
        current_b = b_start_points.copy()
        integrals = []
        d = 0
        while True:
            while current_point[d] + div_widths[d] >= end_points[d] + div_widths[d]/2:
                d += 1
                if d == dimension_number:
                    res = coefficient * (self.func(*current_a) + 4 * self.func(*current_point) + self.func(*current_b))
                    integrals.append(res)
                    total_integral = sum(integrals)
                    return total_integral
            integrals.append(
                coefficient * (self.func(*current_a) + 4 * self.func(*current_point) + self.func(*current_b)))
            current_point[d] += div_widths[d]
            current_a[d] += div_widths[d]
            current_b[d] += div_widths[d]
            for i in range(d):
                current_point[i] = start_points[i]
                current_a[i] = a_start_points[i]
                current_b[i] = b_start_points[i]
            d = 0

    def monte_carlo(self, n=1000):
        """
        Estimates integral using the Monte Carlo method.
        Parameters
        ----------
        n : int
            Number of samples
        Returns
        -------
        Estimate for the integral
        """
        random_func_values = []
        for i in range(0, n):
            # Use random uniform to generate a random x and evaluate the function at that x
            random_x = random.uniform(self.a[0], self.b[0])
            random_func_values.append(self.func(random_x))
        mean = sum(random_func_values)/len(random_func_values)
        integral = mean * (self.b[0] - self.a[0])
        return integral

    def monte_carlo_n_dim(self, n):
        """
        Monte Carlo method generalised to n dimensions.
        Parameters
        ----------
        n : int
            Number of samples
        Returns
        -------
        Estimate for the integral.
        """
        dimension_widths = []
        dimension_number = len(self.a)
        for i in range(dimension_number):
            dimension_width = self.b[i] - self.a[i]
            dimension_widths.append(dimension_width)
        widths_product = np.prod(dimension_widths)
        rand_func_vals = []
        for i in range(n):
            coords = []
            for j in range(dimension_number):
                coords.append(random.uniform(self.a[j], self.b[j]))
            rand_func_vals.append(self.func(*coords))
        mean = sum(rand_func_vals) / n
        return mean * widths_product
