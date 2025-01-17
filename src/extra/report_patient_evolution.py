#!/usr/bin/python3

import sys
import argparse
import numpy as np
import delivery_timing
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# plt.style.use('ggplot')
mpl.rcParams['axes.labelpad'] = 2
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.size'] = 6
mpl.rcParams['figure.titlesize'] = 14
mpl.rcParams['figure.dpi'] = 250
mpl.rcParams['figure.figsize'] = 16, 7


def my_max(l):
    m = -100000
    for i in l:
        b = i.max()
        if b > m:
            m = b
    return m


def my_min(l):
    m = 100000
    for i in l:
        b = i.min()
        if b < m:
            m = b
    return m


def find_range(a, margin=5., va=None):
    if va is not None:
        for i in range(len(a)):
            a[i] += va[i]
    a_min = my_min(a)
    a_max = my_max(a)
    rng = a_max - a_min
    return [a_min - margin/100.*rng, a_max + margin/100.*rng]


def is_outlier(value, p25, p75):
    """Check if value is an outlier
    """
    lower = p25 - 1.5 * (p75 - p25)
    upper = p75 + 1.5 * (p75 - p25)
    return value <= lower or value >= upper
 
 
def get_indices_of_outliers(values):
    """Get outlier indices (if any)
    """
    p25 = np.percentile(values, 25)
    p75 = np.percentile(values, 75)
    IQR = 1.5 * (p75 - p25)
    indices_of_outliers = list()
    for ind, value in enumerate(values):
        if is_outlier(value, p25, p75):
            indices_of_outliers.append(ind)
    return indices_of_outliers, IQR


def add_ellipse(ax, x, y, width, height, colr):
    ax.add_patch(
        patches.Ellipse(
            (x, y), width, height,
            facecolor="None",
            edgecolor=colr, linewidth=0.1
        )
    )


def add_rectangle(ax, x, y, width, height, color, alpha):
    ax.add_patch(
        patches.Rectangle(
            (x, y), width, height,
            facecolor=color, alpha=alpha,
            edgecolor='black', linewidth=0.1
        )
    )


def plot_quiver(ax, x, y, vx, vy, npoints, colors, beamid, spotid, outliers):
    dummy = vx if vx.any() or vy.any() else np.full((npoints, 1), 0.00000001)
    ax.quiver(x, y, dummy, vy, color=colors, angles='xy', scale_units='xy', scale=1, pivot='tail', alpha=1)
    # Plot outliers
    width = 0.05*(x.max()-x.min())
    height = 0.05*(y.max()-y.min())
    for i in outliers:
        add_ellipse(ax, x[i]+0.5*dummy[i], y[i]+0.5*vy[i], width, height, 'red')
        ax.annotate('B' + str(beamid[i]) + 'S' + str(spotid[i]),
                    xy=(x[i], y[i]), fontsize=6,
                    xytext=(x[i], y[i]))


def get_vf_file_data(file):
    r = np.genfromtxt(file, skip_header=1, delimiter=' ',
                      names=['vx', 'vy', 'vz', 'x', 'y', 'z', 'beamid', 'spotid']).T
    vx = np.array(r['vx'])
    vy = np.array(r['vy'])
    vz = np.array(r['vz'])
    x = np.array(r['x'])
    y = np.array(r['y'])
    z = np.array(r['z'])
    beamid = np.array(r['beamid']).astype(int)
    spotid = np.array(r['spotid']).astype(int)

    x = x if hasattr(x, "__len__") else np.array([x])
    y = y if hasattr(y, "__len__") else np.array([y])
    z = z if hasattr(z, "__len__") else np.array([z])
    vx = vx if hasattr(vx, "__len__") else np.array([vx])
    vy = vy if hasattr(vy, "__len__") else np.array([vy])
    vz = vz if hasattr(vz, "__len__") else np.array([vz])
    beamid = beamid if hasattr(beamid, "__len__") else np.array([beamid])
    spotid = spotid if hasattr(spotid, "__len__") else np.array([spotid])

    npoints = len(x)

    d = np.sqrt(vx*vx + vy*vy + vz*vz)
    ang_x = np.array([np.arccos(vx[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])
    ang_y = np.array([np.arccos(vy[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])
    ang_z = np.array([np.arccos(vz[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])

    return vx, vy, vz, x, y, z, beamid, spotid, npoints, d, ang_x, ang_y, ang_z


def analize_vf(vf_files, pp):
    nfracs = len(vf_files)

    vx = list()
    vy = list()
    vz = list()
    x = list()
    y = list()
    z = list()
    beamid = list()
    spotid = list()
    npoints = list()
    d = list()
    ang_x = list()
    ang_y = list()
    ang_z = list()
    for ifrac, vf_file in enumerate(vf_files):
        vx_, vy_, vz_, x_, y_, z_, beamid_, spotid_, npoints_, d_, ang_x_, ang_y_, ang_z_ = get_vf_file_data(vf_file)
        vx.append(vx_)
        vy.append(vy_)
        vz.append(vz_)
        x.append(x_)
        y.append(y_)
        z.append(z_)
        beamid.append(beamid_)
        spotid.append(spotid_)
        npoints.append(npoints_)
        d.append(d_)
        ang_x.append(ang_x_)
        ang_y.append(ang_y_)
        ang_z.append(ang_z_)

    fig = plt.figure()

    fig.suptitle('Vector field size')
    nbins = 50
    color_range = find_range(d, margin=0.)
    xrange = find_range(d)
    ymax_list = list()
    ax_list = list()
    cm = plt.cm.get_cmap('rainbow')
    for ifrac in range(nfracs):
        ax = fig.add_subplot(4, nfracs, ifrac+1)
        if d[ifrac].max() == d[ifrac].min() or xrange[1] == xrange[0]:
            this_nbins = nbins
        else:
            this_nbins = int(1+np.ceil(nbins * (d[ifrac].max() - d[ifrac].min()) / (xrange[1] - xrange[0])))
        bins_y, bins_x = np.histogram(d[ifrac], this_nbins)
        ymax_list.append(bins_y.max())
        hist_colors = [cm((i-color_range[0])/(color_range[1]-color_range[0])) for i in bins_x]
        ax.bar(bins_x[:-1], bins_y, color=hist_colors, width=bins_x[1]-bins_x[0], alpha=0.75, linewidth=0)
        ax.set_xlim(xrange)
        ax.set_xlabel('Vector size (mm)', fontsize=8)
        ax_list.append(ax)
    for ax in ax_list:
        ax.set_ylim([0, 1.05*max(ymax_list)])
    pp.savefig(fig, bbox_inches='tight')
    fig.clf()

    fig.suptitle('Angular histograms')
    nangles = int(360/4)
    rmax_list = list()
    ax_list = list()
    for ifrac in range(nfracs):
        ax = fig.add_subplot(3, nfracs, ifrac+1, projection='polar')
        b, _, _ = ax.hist(ang_x[ifrac], nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
        ax.set_title('Angle x', fontsize=10)
        rmax_list.append(max(b))
        ax_list.append(ax)
        ax = fig.add_subplot(3, nfracs, (ifrac+1)+nfracs, projection='polar')
        b, _, _ = ax.hist(ang_y[ifrac], nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
        ax.set_title('Angle y', fontsize=10)
        rmax_list.append(max(b))
        ax_list.append(ax)
        ax = fig.add_subplot(3, nfracs, (ifrac+1)+2*nfracs, projection='polar')
        b, _, _ = ax.hist(ang_z[ifrac], nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
        ax.set_title('Angle z', fontsize=10)
        rmax_list.append(max(b))
        ax_list.append(ax)
    for ax in ax_list:
        ax.set_rmax(np.round(1.05*max(rmax_list)))
        ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1]) * max(rmax_list)))
    pp.savefig(fig, bbox_inches='tight')
    fig.clf()

    fig.suptitle('Vector field at endpoints')
    range_x = find_range(x, 5., vx)
    range_y = find_range(y, 5., vy)
    range_z = find_range(z, 5., vz)
    color_range = find_range(d, 0.)
    cm = plt.cm.get_cmap('rainbow')
    for ifrac in range(nfracs):
        # Detect outliers
        outliers_x, IQR_x = get_indices_of_outliers(x[ifrac])
        outliers_y, IQR_y = get_indices_of_outliers(y[ifrac])
        outliers_z, IQR_z = get_indices_of_outliers(z[ifrac])
        outliers = np.round(np.unique(np.concatenate((outliers_x, outliers_y, outliers_z)))).astype(int)

        hist_colors = [cm((i - color_range[0]) / (color_range[1] - color_range[0])) for i in d[ifrac]]

        ax = fig.add_subplot(3, nfracs, ifrac+1)
        plot_quiver(ax, x[ifrac], y[ifrac], vx[ifrac], vy[ifrac], npoints[ifrac],
                    hist_colors, beamid[ifrac], spotid[ifrac], outliers)
        ax.set_xlabel('pos x (mm)', fontsize=8)
        ax.set_ylabel('pos y (mm)', fontsize=8)
        ax.set_xlim(range_x)
        ax.set_ylim(range_y)

        ax = fig.add_subplot(3, nfracs, (ifrac+1)+nfracs)
        plot_quiver(ax, y[ifrac], z[ifrac], vy[ifrac], vz[ifrac], npoints[ifrac],
                    hist_colors, beamid[ifrac], spotid[ifrac], outliers)
        ax.set_xlabel('pos y (mm)', fontsize=8)
        ax.set_ylabel('pos z (mm)', fontsize=8)
        ax.set_xlim(range_y)
        ax.set_ylim(range_z)

        ax = fig.add_subplot(3, nfracs, (ifrac+1)+2*nfracs)
        plot_quiver(ax, z[ifrac], x[ifrac], vz[ifrac], vx[ifrac], npoints[ifrac],
                    hist_colors, beamid[ifrac], spotid[ifrac], outliers)
        ax.set_xlabel('pos z (mm)', fontsize=8)
        ax.set_ylabel('pos x (mm)', fontsize=8)
        ax.set_xlim(range_z)
        ax.set_ylim(range_x)

    pp.savefig(fig, bbox_inches='tight')
    fig.clf()


def get_shifts_file_data(shifts_file):
    r = np.genfromtxt(shifts_file, skip_header=1, delimiter=' ',
                      names=['e', 'x', 'y', 'z', 'd', 'beamid', 'spotid']).T
    all_de = np.array(r['e']/1e6)
    all_x = 10*np.array(r['x'])
    all_y = 10*np.array(r['y'])
    # all_z  = np.array([r['z']])
    # all_d  = np.array([r['d']])
    beamid = np.array(r['beamid']).astype(int)
    # spotid = [r['spotid']]
    # Assure it has length
    all_de = all_de if hasattr(all_de, "__len__") else np.array([all_de])
    all_x = all_x if hasattr(all_x, "__len__") else np.array([all_x])
    all_y = all_y if hasattr(all_y, "__len__") else np.array([all_y])
    beamid = beamid if hasattr(beamid, "__len__") else np.array([beamid])
    # Round
    all_de = np.round(all_de, decimals=3)
    all_x = np.round(all_x, decimals=3)
    all_y = np.round(all_y, decimals=3)
    return all_de, all_x, all_y, beamid


def get_tramp_file_data(tramp_file):
    tramp_r = np.genfromtxt(tramp_file, skip_header=12, names=['e', 'x', 'y', 'w']).T
    tramp_e = np.array(tramp_r['e'])
    tramp_x = np.array(tramp_r['x'])
    tramp_y = np.array(tramp_r['y'])
    tramp_w = np.array(tramp_r['w'])
    tramp_e = tramp_e if hasattr(tramp_e, "__len__") else np.array([tramp_e])
    tramp_x = tramp_x if hasattr(tramp_x, "__len__") else np.array([tramp_x])
    tramp_y = tramp_y if hasattr(tramp_y, "__len__") else np.array([tramp_y])
    tramp_w = tramp_w if hasattr(tramp_w, "__len__") else np.array([tramp_w])
    # Round
    tramp_e = np.round(tramp_e, decimals=3)
    tramp_x = np.round(tramp_x, decimals=3)
    tramp_y = np.round(tramp_y, decimals=3)
    return tramp_e, tramp_x, tramp_y, tramp_w


def analize_tramp(shifts_files, tramp_files, pp):
    nfracs = len(shifts_files)
    if not isinstance(tramp_files, (list, tuple)):
        tramp_files = np.array([tramp_files])
    nbeams = len(tramp_files)

    all_de = list()
    all_x = list()
    all_y = list()
    beamid = list()
    for ifrac, shifts_file in enumerate(shifts_files):
        all_de_, all_x_, all_y_, beamid_ = get_shifts_file_data(shifts_file)
        all_de.append(all_de_)
        all_x.append(all_x_)
        all_y.append(all_y_)
        beamid.append(beamid_)

    for ibeam in range(nbeams):
        tramp_file = tramp_files[ibeam]
        tramp_e, tramp_x, tramp_y, tramp_w = get_tramp_file_data(tramp_file)
        unique_energies = np.unique(tramp_e)[::-1]
        number_layers = len(unique_energies)
        layer_id = np.array([int(np.squeeze(np.where(unique_energies == i))) for i in tramp_e])

        de = list()
        x = list()
        y = list()
        vx = list()
        vy = list()
        d = list()
        for ifrac in range(nfracs):
            de.append(np.array(all_de[ifrac][beamid[ifrac] == ibeam]))
            x.append(np.array(all_x[ifrac][beamid[ifrac] == ibeam]))
            y.append(np.array(all_y[ifrac][beamid[ifrac] == ibeam]))
            vx.append(x[ifrac] - tramp_x)
            vy.append(y[ifrac] - tramp_y)
            d.append(np.sqrt(vx[ifrac]*vx[ifrac] + vy[ifrac]*vy[ifrac]))

        fig = plt.figure()

        # FIGURE 1, PLOT 1 --------------------------------
        fig.suptitle('Beam: {}. Energy layer distribution.'.format(ibeam))
        for ifrac in range(nfracs):
            ax = fig.add_subplot(3, nfracs, ifrac+1)
            accu_len = 0
            for i, layer_energy in enumerate(unique_energies):
                bool_mask = tramp_e == layer_energy
                layer = (tramp_e+de[ifrac])[bool_mask]
                average = np.mean(layer)
                pos_text_y = 1.01*max(layer) if max(layer) > 1.05*layer_energy else 1.01*1.05*layer_energy
                ax.annotate("{:.2f} %".format(
                            100*(average-layer_energy)/layer_energy), xy=(accu_len, pos_text_y), fontsize=4)
                box_x0 = accu_len-1 if i != 0 else accu_len
                box_x1 = len(layer) if i != 0 else len(layer)-1
                add_rectangle(ax, box_x0, min(layer), box_x1, max(layer)-min(layer), 'blue', 0.1)
                add_rectangle(ax, box_x0, 0.95*layer_energy, box_x1, 0.1*layer_energy, 'green', 0.075)
                accu_len += len(layer)
            ax.step(range(len(tramp_e)), tramp_e, color='black', alpha=0.9)
            ax.plot(np.sort(tramp_e+de[ifrac])[::-1], linestyle='None', marker='o', color='black',
                    markersize=1, alpha=0.4, markeredgecolor='black', markeredgewidth=0.1)
            ax.plot(tramp_e+de[ifrac], linestyle='None', color='red',
                    marker='o', markersize=1.75, markeredgewidth=0.25)
            ax.set_xlabel('Spot number', fontsize=7)
            ax.set_ylabel('Energy (MeV)', fontsize=7)
            ax.annotate('N. spots = ' + str(len(tramp_e)) + '\nOrig. layers  = ' +
                        str(number_layers) + '\nAdapt. layers = ' +
                        str(len(np.unique(tramp_e+de[ifrac]))),
                        xy=(5, min(min(tramp_e+de[ifrac]), min(tramp_e))), fontsize=6,
                        textcoords='axes fraction', xytext=(0.04, 0.04))
        pp.savefig(fig, bbox_inches='tight')
        fig.clf()

        # FIGURE 1, PLOT 2 --------------------------------
        fig.suptitle('Beam: {}. Energy layer shifts.'.format(ibeam))
        nbins = 50
        color_range = [my_min(de), my_max(de)]
        xrange = find_range(de)
        ymax_list = list()
        ax_list = list()
        cm = plt.cm.get_cmap('rainbow')
        for ifrac in range(nfracs):
            ax = fig.add_subplot(3, nfracs, ifrac + 1)
            this_nbins = int(np.ceil(nbins*(de[ifrac].max() - de[ifrac].min())/(xrange[1] - xrange[0])))
            bins_y, bins_x = np.histogram(de[ifrac], this_nbins)
            ymax_list.append(bins_y.max())
            hist_colors = [cm((i - color_range[0]) / (color_range[1] - color_range[0])) for i in bins_x]
            ax.bar(bins_x[:-1], bins_y, color=hist_colors, width=bins_x[1] - bins_x[0], alpha=0.75, linewidth=0)
            ax.set_xlim(xrange)
            ax.set_xlabel('Energy shift (MeV)', fontsize=8)
            ax_list.append(ax)
        for ax in ax_list:
            ax.set_ylim([0, 1.05 * max(ymax_list)])
        pp.savefig(fig, bbox_inches='tight')
        fig.clf()

        # FIGURE 1, PLOT 3 --------------------------------
        fig.suptitle('Beam: {}. Source plane positions.'.format(ibeam))
        cmap = plt.cm.get_cmap('rainbow')
        color_range = find_range(de, 0.)
        for ifrac in range(nfracs):
            ax = fig.add_subplot(3, nfracs, ifrac+1)
            for i in range(len(tramp_x)):
                if d[ifrac][i] > 0.001:
                    ax.arrow(0.1 * vx[ifrac][i] + tramp_x[i], 0.1 * vy[ifrac][i] + tramp_y[i],
                             0.8 * vx[ifrac][i], 0.8 * vy[ifrac][i],
                             fc='k', ec='k', alpha=0.25, zorder=0)
            scat_colors = [cm((i - color_range[0]) / (color_range[1] - color_range[0])) for i in de[ifrac]]
            # scat_colors = np.array(cmap(layer_id / max(layer_id.max(), len(layer_id))))
            ax.scatter(tramp_x, tramp_y, s=10, linewidth=0.5, zorder=1,
                       edgecolors='black', alpha=0.75, facecolors='')
            ax.scatter(x[ifrac], y[ifrac], s=10, linewidth=0.5, zorder=2, edgecolors='k',
                       facecolors=scat_colors)
            ax.set_xlabel('X (mm)', fontsize=7)
            ax.set_ylabel('Y (mm)', fontsize=7)
        pp.savefig(fig, bbox_inches='tight')
        fig.clf()

        # FIGURE 1, PLOT 4 --------------------------------
        fig.suptitle('Beam: {}. Slow direction layer shifts.'.format(ibeam))
        nbins = 50
        color_range = find_range(d, 0.)
        xrange = find_range(d)
        ymax_list = list()
        ax_list = list()
        cm = plt.cm.get_cmap('rainbow')
        for ifrac in range(nfracs):
            ax = fig.add_subplot(3, nfracs, ifrac + 1)
            this_nbins = int(np.ceil(nbins*(d[ifrac].max() - d[ifrac].min())/(xrange[1] - xrange[0])))
            bins_y, bins_x = np.histogram(d[ifrac], this_nbins)
            ymax_list.append(bins_y.max())
            hist_colors = [cm((i - color_range[0]) / (color_range[1] - color_range[0])) for i in bins_x]
            ax.bar(bins_x[:-1], bins_y, color=hist_colors, width=bins_x[1] - bins_x[0], alpha=0.75, linewidth=0)
            ax.set_xlim(xrange)
            ax.set_xlabel('Position shift (mm)', fontsize=8)
            ax_list.append(ax)
        for ax in ax_list:
            ax.set_ylim([0, 1.05 * max(ymax_list)])
        pp.savefig(fig, bbox_inches='tight')
        fig.clf()

        # FIGURE 1, PLOT 5 --------------------------------
        fig.suptitle('Beam: {}. Slow direction layer positions.'.format(ibeam))
        for ifrac in range(nfracs):
            ax = fig.add_subplot(3, nfracs, ifrac+1)
            ax.step(range(len(tramp_y)), tramp_y, color='black', alpha=0.5)
            ax.step(range(len(tramp_y)), y[ifrac], color='blue', alpha=0.5)
            ax.plot(y[ifrac], linestyle='None', linewidth=0.5, color='red', marker='o', markersize=2)
            ax.set_xlabel('Spot number', fontsize=7)
            ax.set_ylabel('Slow dimension pos (mm)', fontsize=7)
            ax.annotate('N. spots = ' + str(len(tramp_e)) + '\n' +
                        'Orig. layers = ' + str(len(np.unique(tramp_y))) + '\n' +
                        'Adapt. layers = ' + str(len(np.unique(y[ifrac]))),
                        xy=(5, min(min(y[ifrac]), min(tramp_y))), fontsize=6,
                        textcoords='axes fraction', xytext=(0.04, 0.04))
        pp.savefig(fig, bbox_inches='tight')
        fig.clf()

        # FIGURE 1, PLOT 6 --------------------------------
        fig.suptitle('Beam: {}. Slow direction layer positions.'.format(ibeam))
        for ifrac in range(nfracs):
            time, summary = delivery_timing.get_timing(tramp_e, tramp_x, tramp_y, tramp_w)
            time_adapted, summary_adapted = delivery_timing.get_timing(tramp_e+de[ifrac], x[ifrac], y[ifrac],
                                                                       tramp_w)
            max_time = 120. if 120. > time[-1] else 1.5*time[-1]  # s
            default_energy_switch = 2.5                           # s
            time_ideal, summary_ideal = delivery_timing.get_timing(tramp_e, x[ifrac], y[ifrac],
                                                                   tramp_w, energy_switch_time=False)
            total_ideal = time_ideal[-1]
            max_layers = np.floor((max_time - total_ideal)/default_energy_switch)
            # print('Maximum theoretical energy layers: {}'.format(max_layers))

            summary = summary.replace('Summary', 'Orig.')
            summary_adapted = summary_adapted.replace('Summary', 'Adapt.')

            ax = fig.add_subplot(3, nfracs, ifrac+1)
            ax.plot(time, color='blue', alpha=0.5)
            ax.plot(time_adapted, color='red')
            accu_len = 0
            for layer_energy in unique_energies:
                time_lyr = time[tramp_e == layer_energy]
                time_adapted_lyr = time_adapted[tramp_e == layer_energy]
                add_rectangle(ax, accu_len, min(time_lyr), len(time_lyr),
                              max(time_adapted_lyr)-min(time_lyr), 'green', 0.1)
                accu_len += len(time_lyr)
            ax.set_xlabel('Spot number', fontsize=7)
            ax.set_ylabel('Time (s)', fontsize=7)
            ax.tick_params(labelsize=6)
            bbox_props = dict(boxstyle="Round", facecolor="blue", alpha=0.1)
            bbox_props_adapt = dict(boxstyle="Round", facecolor="red", alpha=0.1)
            ax.annotate(summary, family='monospace',
                        textcoords='axes fraction', xytext=(0.03, 0.64),
                        xy=(5, min(time)), fontsize=4, bbox=bbox_props)
            ax.annotate(summary_adapted, family='monospace',
                        textcoords='axes fraction', xytext=(0.31, 0.64),
                        xy=(5, min(time)), fontsize=4, bbox=bbox_props_adapt)
        pp.savefig(fig, bbox_inches='tight')
        fig.clf()


def main(argv):
    # Define the top-level parser
    parser = argparse.ArgumentParser(description='Analysis of adaptation')
    # Parse module:
    parser.add_argument('--vf',     nargs='+', help='File with vector field values and locations', required=True)
    parser.add_argument('--shifts', nargs='+', help='File with energy and locations shifts', required=True)
    parser.add_argument('--tramps', nargs='+', help='Tramp files of the original plan', required=True)
    parser.add_argument('--outdir', help='Directory to output analysis', default='./')
    parser.add_argument('--outfile', help='Report file name. It will be forced to be a pdf.',
                        default='adapt_report.pdf')

    args = parser.parse_args(argv)

    if not args.outfile.endswith('.pdf'):
        args.outfile += '.pdf'
    outfile = args.outdir + "/" + args.outfile
    pp = PdfPages(outfile)
    analize_vf(args.vf, pp)
    analize_tramp(args.shifts, args.tramps, pp)
    pp.close()


if __name__ == "__main__":
    main(sys.argv[1:])
