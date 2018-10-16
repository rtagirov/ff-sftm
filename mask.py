import numpy as np

def spot_patch(x, y, x_min, x_max, y_min, y_max, grid_max):

    yc = (y >= y_min) & (y <= y_max)

    if x_min >= 0.0 and x_max <= grid_max:

        xc = (x >= x_min) & (x <= x_max)

        return x[np.where(xc & yc)], y[np.where(xc & yc)]

    if x_min < 0.0:

        if x_max > 0.0:

            x_min = grid_max + x_min

            xc = ((x <= x_max) | (x >= x_min))

            return x[np.where(xc & yc)], y[np.where(xc & yc)]

        if x_max < 0.0:

            x_min = grid_max + x_min
            x_max = grid_max + x_max

            xc = (x >= x_min) & (x <= x_max)

            return x[np.where(xc & yc)], y[np.where(xc & yc)]

    if x_max > grid_max:

        if x_min < grid_max:

            x_max = x_max - grid_max

            xc = ((x <= x_max) | (x >= x_min))

            return x[np.where(xc & yc)], y[np.where(xc & yc)]

        if x_min > grid_max:

            x_min = x_min - grid_max
            x_max = x_max - grid_max

            xc = (x >= x_min) & (x <= x_max)

            return x[np.where(xc & yc)], y[np.where(xc & yc)]

