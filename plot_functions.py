import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def multiplot(data=None, scan_range=None, labels=None, title_list=None, scale='linear', color_map='viridis', interpolation='spline36'):
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
    print('rows ', rows, ' cols ', cols)
    axes = []
    fig = plt.figure()
    for k in range(num_plots):
        axes.append(fig.add_subplot(rows, cols, k + 1))
        subplot_title = (title_list[k])
        axes[-1].set_title(subplot_title)
        if scale == 'log':
            im = plt.imshow(data[k], cmap=color_map, norm=colors.LogNorm(), origin='lower', interpolation=interpolation,
                       extent=scan_range)
        else:
            im = plt.imshow(data[k], cmap=color_map, origin='lower', interpolation=interpolation, extent=scan_range)
        if labels:
            plt.xlabel(labels[0])
            plt.ylabel(labels[1])
        plt.colorbar(im)

    fig.tight_layout()
    plt.show()

    return


def plot(data=None, scan_range=None, labels=None, title=None, scale='linear', color_map='viridis', interpolation='spline36'):
    if data is None:
        print('Nothing to plot, kindly provide the data')
        return
    if scan_range is None:
        print('Scan range not given')
        scan_range = [0, 1, 0, 1]

    plt.figure()
    if scale == 'log':
        plt.imshow(data, cmap=color_map, norm=colors.LogNorm(), origin='lower', interpolation=interpolation,
                       extent=scan_range)
    else:
        plt.imshow(data, cmap=color_map, origin='lower', interpolation=interpolation, extent=scan_range)
    plt.colorbar()
    if title:
        plt.title(title)
    if labels:
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
    plt.show()

    return

