import numpy as np
import matplotlib.pyplot as plt


def multiplot(data=None, scan_range=None, labels=None, title_list=None, scale='linear', color_map='viridis', interpolation='spline36'):
    """
    Plot multiple dataset for spectral and evolution data
    :param data: List of spectra or dipole expectation values or any other variable of interest
    :param scan_range: The min and max of both axis in the format [xmin, xmax, ymin, ymax]
    :param labels: List of label for each axis
    :param title_list: List of titles for each plot
    :param scale: Scaling of the data points, two choices are 'linear' and 'log'
    :param color_map: Choice of colormap
    :param interpolation: Interpolation for points in plot.
    :return: Does not return anything
    """
    if data is None:
        print('Nothing to plot, kindly provide the data')
        return
    if scan_range is None:
        print('Scan range not given')
        scan_range = [0, 1, 0, 1]
    if title_list is None:
        print('titles not given')
        title_list = [str(x + 1) for x in range(len(data))]

    num_plots = len(data)  # number of plots (depends on the length of data list)
    if num_plots <= 3:
        rows = 1
        cols = num_plots
    else:
        rows = int(np.ceil(num_plots / 3))
        cols = 3

    if scale == 'log':
        data = np.array([log_scale(s.real) for s in data])

    axes = []
    fig = plt.figure(figsize=(16, 4))
    for k in range(num_plots):
        axes.append(fig.add_subplot(rows, cols, k + 1))
        subplot_title = (title_list[k])
        axes[-1].set_title(subplot_title)
        im = plt.imshow(data[k], cmap=color_map, origin='lower', interpolation=interpolation, extent=scan_range, aspect=1)
        if labels:
            plt.xlabel(labels[0])
            plt.ylabel(labels[1])
        plt.colorbar(im, ax=axes[-1])

    fig.tight_layout()
    plt.show()

    return


def plot(data=None, scan_range=None, labels=None, title=None, scale='linear', color_map='viridis', interpolation='spline36'):
    """
    Plot multiple dataset for spectral and evolution data
    :param data: List of spectra or dipole expectation values or any other variable of interest
    :param scan_range: The min and max of both axis in the format [xmin, xmax, ymin, ymax]
    :param labels: List of label for each axis
    :param title: Title the plot
    :param scale: Scaling of the data points, two choices are 'linear' and 'log'
    :param color_map: Choice of colormap
    :param interpolation: Interpolation for points in plot.
    :return: Does not return anything
    """
    if data is None:
        print('Nothing to plot, kindly provide the data')
        return
    if scan_range is None:
        print('Scan range not given')
        scan_range = [0, 1, 0, 1]

    plt.figure()
    if scale == 'log':
        data = log_scale(data.real)
    plt.imshow(data, cmap=color_map, origin='lower', interpolation=interpolation, extent=scan_range, aspect='auto')
    plt.colorbar()
    if title:
        plt.title(title)
    if labels:
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
    plt.show()

    return


def log_scale(z):
    """
    Simple function for rescaling the 2D input matrix to log scale.
    Warning: the numbers between -1 and 1 are mapped to zero.
    """
    x, y = np.shape(z)
    for n in range(x):
        for m in range(y):
            if z[n, m] >= 1:
                z[n, m] = np.log(z[n, m])
            elif z[n, m] <= -1:
                z[n, m] = -np.log(-z[n, m])
            else:
                z[n, m] = 0
    return z
