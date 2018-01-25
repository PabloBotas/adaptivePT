#!/usr/bin/python3

import sys
import argparse
import numpy as np
import delivery_timing
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
# plt.style.use('ggplot')
mpl.rcParams['axes.labelpad'] = 2
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.size'] = 6
mpl.rcParams['figure.titlesize'] = 14
mpl.rcParams['figure.dpi'] = 250
mpl.rc('text', usetex=True)

def find_range(a, margin=5., va=None):
    if va is not None:
        a += va
    a_min = a.min()
    a_max = a.max()
    rng = a_max - a_min
    return [a_min - margin/100.*rng, a_max + margin/100.*rng]


def is_outlier(value, p0, p1):
    """Check if value is an outlier
    """
    IQR = (p1 - p0)
    lower = p0 - 1.5*IQR
    upper = p1 + 1.5*IQR
    return value <= lower or value >= upper
 
 
def get_indices_of_outliers(values):
    """Get outlier indices (if any)
    """
    p0 = np.percentile(values, 10)
    p1 = np.percentile(values, 90)
     
    indices_of_outliers = []
    for ind, value in enumerate(values):
        if is_outlier(value, p0, p1):
            indices_of_outliers.append(ind)
    return indices_of_outliers


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

def add_axis_frame(fig, ax, color, alpha=1):
    ax.relim()
    x = ax.get_xbound()
    y = ax.get_ybound()
    print(x,y)
    rect = patches.Rectangle((x[0], 5.28),
                             x[1]-x[0], 12.8-5.28,
                             edgecolor='red', alpha=0.9,
                             linewidth=4,
                             fill=False)
    rect.set_clip_on(False)
    ax.add_patch(rect)

def freedman_diaconis_bins(x):
    r = find_range(x)
    perc = np.percentile(x, [75, 25])
    if r[0] == r[1] or perc[0] == perc[1]:
        return 1
    return int((r[1] - r[0])/(2*np.subtract(*perc)/np.power(len(x), 1/3)) + 1)


def shimazaki_bins(x):
    from scipy.optimize import minimize
    r = np.max(x)-np.min(x)
    delta = r/10
    bnds = ((r/1000, None),)
    def cost(delta):
        # print(delta[0], int(r/delta[0] + 1), r)
        k,_ = np.histogram(x, int(r/delta[0] + 1))
        ave = np.mean(k)
        v = np.var(k)
        return (2*ave-v)/(delta[0]*delta[0])
    result = minimize(cost, ((delta),), tol=1e-6, bounds=bnds)
    return int(r/result.x+1)


def optimal_n_bins(x):
    return freedman_diaconis_bins(x)


def analize_vf(vf_file, pp):
    r = np.genfromtxt(vf_file, skip_header=1, delimiter=' ',
                      names=['vx', 'vy', 'vz', 'x', 'y', 'z', 'beamid', 'spotid', 'gantryangle', 'couchangle']).T
    vx = np.array(r['vx'])
    vy = np.array(r['vy'])
    vz = np.array(r['vz'])
    x = np.array(r['x'])
    y = np.array(r['y'])
    z = np.array(r['z'])
    beamid = np.array(r['beamid']).astype(int)
    spotid = np.array(r['spotid']).astype(int)
    beamangle = np.array(r['gantryangle'])
    couchangle = np.array(r['couchangle'])

    x = x if x.size > 1 else np.array([x])
    y = y if y.size > 1 else np.array([y])
    z = z if z.size > 1 else np.array([z])
    vx = vx if vx.size > 1 else np.array([vx])
    vy = vy if vy.size > 1 else np.array([vy])
    vz = vz if vz.size > 1 else np.array([vz])
    beamid = beamid if beamid.size > 1 else np.array([beamid])
    spotid = spotid if spotid.size > 1 else np.array([spotid])
    beamangle = beamangle if beamangle.size > 1 else np.array([beamangle])
    couchangle = couchangle if couchangle.size > 1 else np.array([couchangle])
    beamangle = np.unique(beamangle)

    npoints = x.size

    d = np.sqrt(vx*vx + vy*vy + vz*vz)
    ang_x = np.array([np.arccos(vx[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])
    ang_y = np.array([np.arccos(vy[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])
    ang_z = np.array([np.arccos(vz[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])

    nbins = optimal_n_bins(d)
    nangles1 = optimal_n_bins(ang_x)
    nangles2 = optimal_n_bins(ang_y)
    nangles3 = optimal_n_bins(ang_z)

    fig = plt.figure(figsize=(10, 6))
    ax1 = plt.subplot2grid((2, 12), (0, 0), colspan=3)
    ax2 = plt.subplot2grid((2, 12), (0, 3), colspan=3, projection='polar')
    ax3 = plt.subplot2grid((2, 12), (0, 6), colspan=3, projection='polar')
    ax4 = plt.subplot2grid((2, 12), (0, 9), colspan=3, projection='polar')
    ax5 = plt.subplot2grid((2, 12), (1, 0), colspan=4)
    ax6 = plt.subplot2grid((2, 12), (1, 4), colspan=4)
    ax7 = plt.subplot2grid((2, 12), (1, 8), colspan=4)
    fig.suptitle('Vector field analysis at endpoints')

    # FIGURE 1, PLOT 1 --------------------------------
    bins_y, bins_x = np.histogram(d, bins=nbins)
    cm = plt.cm.get_cmap('rainbow')
    hist_colors = [cm(((i-bins_x.min())/(bins_x.max()-bins_x.min()))) for i in bins_x]
    ax1.bar(bins_x[:-1], bins_y, color=hist_colors, width=bins_x[1]-bins_x[0], alpha=0.75)
    ax1.set_xlabel('Vector size (cm)', fontsize=8)

    # FIGURE 1, POLAR PLOT 1 --------------------------------
    weights = d/float(np.sum(d)) if np.any(d) else np.ones_like(d)
    b, _, _ = ax2.hist(ang_x, nangles1, weights=weights, histtype='step',
                       alpha=1, color='r')
    ax2.set_rorigin(-0.4*max(b))
    ax2.tick_params(direction='out', pad=-2)
    ax2.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b), decimals=2))
    ax2.set_title('Angle x', fontsize=11)

    # FIGURE 1, POLAR PLOT 2 --------------------------------
    b, _, _ = ax3.hist(ang_y, nangles2, weights=weights, histtype='step',
                      alpha=1, color='r')
    color_list = ['b', 'g', 'c', 'm', 'y', 'k']
    for i,ang in enumerate(beamangle):
        temp = ang + 270.0*np.pi/180
        print((temp*180.0/np.pi) % 360.0)
        ax3.axvline(x=temp, linewidth=20, alpha=0.2, label='beam '+str(i),
                   zorder=0, color=color_list[i])
    ax3.set_rorigin(-0.4*max(b))
    ax3.tick_params(direction='out', pad=-2)
    ax3.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b), decimals=2))
    ax3.set_title('Angle y', fontsize=11)
    leg = ax3.legend(shadow=True, fancybox=True)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)

    # FIGURE 1, POLAR PLOT 3 --------------------------------
    b, _, _ = ax4.hist(ang_z, nangles3, weights=weights, histtype='step',
                      alpha=1, color='r')
    ax4.set_rorigin(-0.4*max(b))
    ax4.tick_params(direction='out', pad=-2)
    ax4.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b), decimals=2))
    ax4.set_title('Angle z', fontsize=11)

    # FIGURE 1, QUIVER PLOT 1 --------------------------------
    # Detect outliers
    outliers_x = get_indices_of_outliers(x)
    outliers_y = get_indices_of_outliers(y)
    outliers_z = get_indices_of_outliers(z)
    outliers = np.round(np.unique(np.concatenate((outliers_x, outliers_y, outliers_z)))).astype(int)

    dummy = vx if vx.any() or vy.any() else np.full((npoints, 1), 0.00000001)
    ax5.quiver(x, y, dummy, vy, d, angles='xy', scale_units='xy', scale=1,
              cmap=plt.cm.get_cmap('rainbow'), pivot='tail', alpha=1)
    ax5.set_aspect('equal', 'datalim')
    ax5.set_xlabel('pos x (cm)', fontsize=8)
    ax5.set_ylabel('pos y (cm)', fontsize=8)
    # Detect outliers
    width = 0.05*(x.max()-x.min())
    height = 0.05*(y.max()-y.min())
    for circ in outliers:
        add_ellipse(ax5, x[circ]+0.5*dummy[circ], y[circ]+0.5*vy[circ], width, height, 'red')
        ax5.annotate('B' + str(beamid[circ]) + 'S' + str(spotid[circ]),
                    xy=(x[circ], y[circ]), fontsize=6,
                    xytext=(x[circ], y[circ]))

    # FIGURE 1, QUIVER PLOT 2 --------------------------------
    dummy = vy if vy.any() or vz.any() else np.full((npoints, 1), 0.00000001)
    ax6.quiver(y, z, dummy, vz, d, angles='xy', scale_units='xy', scale=1,
              cmap=plt.cm.get_cmap('rainbow'), pivot='tail', alpha=1)
    ax6.set_aspect('equal', 'datalim')
    ax6.set_xlabel('pos y (cm)', fontsize=8)
    ax6.set_ylabel('pos z (cm)', fontsize=8)
    # Detect outliers
    height = 0.05*(z.max()-z.min())
    width = 0.05*(y.max()-y.min())
    for circ in outliers:
        add_ellipse(ax6, y[circ]+0.5*dummy[circ], z[circ]+0.5*vz[circ], width, height, 'red')
        ax6.annotate('B' + str(beamid[circ]) + 'S' + str(spotid[circ]),
                    xy=(y[circ], z[circ]), fontsize=6,
                    xytext=(y[circ], z[circ]))

    # FIGURE 1, QUIVER PLOT 3 --------------------------------
    dummy = vz if vz.any() or vx.any() else np.full((npoints, 1), 0.00000001)
    ax7.quiver(z, x, dummy, vx, d, angles='xy', scale_units='xy', scale=1,
              cmap=plt.cm.get_cmap('rainbow'), pivot='tail', alpha=1)
    ax7.set_aspect('equal', 'datalim')
    ax7.set_xlabel('pos z (cm)', fontsize=8)
    ax7.set_ylabel('pos x (cm)', fontsize=8)
    # Detect outliers
    width = 0.05*(z.max()-z.min())
    height = 0.05*(x.max()-x.min())
    for circ in outliers:
        add_ellipse(ax7, z[circ]+0.5*dummy[circ], x[circ]+0.5*vx[circ], width, height, 'red')
        ax7.annotate('B' + str(beamid[circ]) + 'S' + str(spotid[circ]),
                    xy=(z[circ], x[circ]), fontsize=6,
                    xytext=(z[circ], x[circ]))

    fig.tight_layout()
    # add_axis_frame(fig, ax5, 'red')
    pp.savefig(fig, bbox_inches='tight', dpi='figure')


def analize_tramp(shifts_file, tramp_files, spots_layer, pp):
    r = np.genfromtxt(shifts_file, comments='#', delimiter=' ',
                      names=['w', 'e', 'x', 'y', 'z', 'd', 'beamid', 'spotid']).T
    all_dw = np.array(r['w'])
    all_de = np.array(r['e']/1e6)
    all_x = 10*np.array(r['x'])
    all_y = 10*np.array(r['y'])
    # all_z  = np.array([r['z']])
    # all_d  = np.array([r['d']])
    beamid = np.array(r['beamid']).astype(int)
    # spotid = [r['spotid']]

    all_dw = all_dw if all_dw.size > 1 else np.array([all_dw])
    all_de = all_de if all_de.size > 1 else np.array([all_de])
    all_x  = all_x  if all_x.size  > 1 else np.array([all_x])
    all_y  = all_y  if all_y.size  > 1 else np.array([all_y])
    beamid = beamid if beamid.size > 1 else np.array([beamid])

    # Round
    all_dw = np.round(all_dw, decimals=8)
    all_de = np.round(all_de, decimals=3)
    all_x = np.round(all_x, decimals=3)
    all_y = np.round(all_y, decimals=3)

    for tramp_num, tramp_file in enumerate(tramp_files, start=0):
        tramp_r = np.genfromtxt(tramp_file, skip_header=12, names=['e', 'x', 'y', 'w']).T
        tramp_e = np.array(tramp_r['e'])
        tramp_x = np.array(tramp_r['x'])
        tramp_y = np.array(tramp_r['y'])
        tramp_w = np.array(tramp_r['w'])

        tramp_e = tramp_e if tramp_e.size > 1 else np.array([tramp_e])
        tramp_x = tramp_x if tramp_x.size > 1 else np.array([tramp_x])
        tramp_y = tramp_y if tramp_y.size > 1 else np.array([tramp_y])
        tramp_w = tramp_w if tramp_w.size > 1 else np.array([tramp_w])

        # Round
        tramp_e = np.round(tramp_e, decimals=3)
        tramp_x = np.round(tramp_x, decimals=3)
        tramp_y = np.round(tramp_y, decimals=3)
        tramp_w = np.round(tramp_w, decimals=8)

        dw = np.array(all_dw[beamid == tramp_num])
        de = np.array(all_de[beamid == tramp_num])
        x = np.array(all_x[beamid == tramp_num])
        y = np.array(all_y[beamid == tramp_num])
        vx = x - tramp_x
        vy = y - tramp_y
        d = np.sqrt(vx*vx + vy*vy)

        unique_energies = np.unique(tramp_e)[::-1]
        number_layers = len(unique_energies)
        layer_id = np.array([int(np.squeeze(np.where(unique_energies == i))) for i in tramp_e])

        fig = plt.figure(figsize = (10, 10))
        fig.suptitle('Tramp adaptations analysis. Beam: {}'.format(tramp_num))

        # FIGURE 1, PLOT 1 --------------------------------
        ax = fig.add_subplot(4, 2, 1)
        accu_len = 0
        for i, layer_energy in enumerate(unique_energies):
            bool_mask = tramp_e == layer_energy
            layer = (tramp_e+de)[bool_mask]
            average = np.mean(layer)
            pos_text_y = 1.01*max(layer) if max(layer) > 1.05*layer_energy else 1.01*1.05*layer_energy
            ax.annotate("{:.2f}".format(
                100*(average-layer_energy)/layer_energy), xy=(accu_len, pos_text_y), fontsize=5)
            box_x0 = accu_len-1 if i != 0 else accu_len
            box_x1 = len(layer) if i != 0 else len(layer)-1
            add_rectangle(ax, box_x0, min(layer), box_x1, max(layer)-min(layer), 'blue', 0.1)
            add_rectangle(ax, box_x0, 0.95*layer_energy, box_x1, 0.1*layer_energy, 'green', 0.075)
            accu_len += len(layer)
        ax.step(range(len(tramp_e)), tramp_e, color='black', alpha=0.9)
        ax.plot(np.sort(tramp_e+de)[::-1], linestyle='None', marker='o', color='black',
                markersize=1, alpha=0.4, markeredgecolor='black', markeredgewidth=0.1)
        ax.plot(tramp_e+de, linestyle='None', color='red', marker='o', alpha=0.75,
                markersize=1.75, markeredgewidth=0.25)
        ax.set_xlabel('Spot number', fontsize=7)
        ax.set_ylabel('Energy (MeV)', fontsize=7)
        ax.annotate('N. spots = ' + str(len(tramp_e)) + '\nOrig. layers  = ' +
                    str(number_layers) + '\nAdapt. layers = ' + str(len(np.unique(tramp_e+de))),
                    xy=(5, min(min(tramp_e+de), min(tramp_e))), fontsize=7,
                    textcoords='axes fraction', xytext=(0.04, 0.04))

        # FIGURE 1, PLOT 2 --------------------------------
        ax = fig.add_subplot(4, 2, 2)
        color_range = find_range(de, 0.)
        if np.abs(color_range[0]) > np.abs(color_range[1]):
            color_range = (-np.abs(color_range[0]), np.abs(color_range[0]))
        else:
            color_range = (-np.abs(color_range[1]), np.abs(color_range[1]))
        if color_range[0] == color_range[1]:
            color_range[0] -= 2
            color_range[0] += 2
        cm = plt.cm.get_cmap('rainbow')
        bins_y, bins_x = np.histogram(de, bins=optimal_n_bins(de))
        hist_colors = [cm((i - color_range[0]) / (color_range[1] - color_range[0])) for i in bins_x]
        ax.bar(bins_x[:-1], bins_y, color=hist_colors, width=bins_x[1] - bins_x[0], alpha=0.75, linewidth=0)
        ax.axvline(x=0, color='k', linewidth=1)
        ax.set_xlim(find_range(de))
        ax.set_xlabel('Energy shift (MeV)', fontsize=7)

        # FIGURE 1, PLOT 3 --------------------------------
        ax = fig.add_subplot(4, 2, 3)
        cmap = plt.cm.get_cmap('rainbow')
        color_range = find_range(de, 0.)
        if np.abs(color_range[0]) > np.abs(color_range[1]):
            color_range = (-np.abs(color_range[0]), np.abs(color_range[0]))
        else:
            color_range = (-np.abs(color_range[1]), np.abs(color_range[1]))
        if color_range[0] == color_range[1]:
            color_range[0] -= 2
            color_range[0] += 2
        for i in range(len(tramp_x)):
            if d[i] > 0.001:
                ax.arrow(0.1 * vx[i] + tramp_x[i], 0.1 * vy[i] + tramp_y[i],
                         0.8 * vx[i], 0.8 * vy[i],
                         fc='k', ec='k', alpha=0.25, zorder=0)
        scat_colors = [cm((i - color_range[0]) / (color_range[1] - color_range[0])) for i in de]
        ax.scatter(tramp_x, tramp_y, s=10, linewidth=0.25, zorder=1,
                   edgecolors='black', alpha=0.75, facecolors='')
        ax.scatter(x, y, s=10, linewidth=0.5, alpha=1, zorder=2, edgecolors='k',
                   facecolors=scat_colors)
        ax.set_xlabel('X (mm)', fontsize=7)
        ax.set_ylabel('Y (mm)', fontsize=7)

        # FIGURE 1, PLOT 4 --------------------------------
        ax = fig.add_subplot(4, 2, 4)
        if d.max()-d.min() < 0.1:
            ax.hist(d, bins=optimal_n_bins(d), range=(d.min()-0.05, d.max()+0.05))
        else:
            ax.hist(d, bins=optimal_n_bins(d))
        ax.set_xlabel('Shift size (mm)', fontsize=7)

        # FIGURE 1, PLOT 5 --------------------------------
        ax = fig.add_subplot(4, 2, 5)
        ax.step(range(len(tramp_y)), tramp_y, color='black', alpha=0.5)
        ax.step(range(len(tramp_y)), y, color='blue', alpha=0.5)
        ax.plot(y, linestyle='None', linewidth=1, color='red', marker='o', markersize=2)
        ax.set_xlabel('Spot number', fontsize=7)
        ax.set_ylabel('Slow dimension pos (mm)', fontsize=7)
        ax.annotate('N. spots = ' + str(len(tramp_e)) + '\n' +
                    'Orig. layers = ' + str(len(np.unique(tramp_y))) + '\n' +
                    'Adapt. layers = ' + str(len(np.unique(y))),
                    xy=(5, min(min(y), min(tramp_y))), fontsize=6,
                    textcoords='axes fraction', xytext=(0.04, 0.04))

        # FIGURE 1, PLOT 6 --------------------------------
        time, summary = delivery_timing.get_timing(tramp_e, tramp_x, tramp_y, tramp_w)
        time_adapted, summary_adapted = delivery_timing.get_timing(tramp_e+de, x, y, tramp_w)
        max_time = 120. if 120. > time[-1] else 1.5*time[-1]  # s
        default_energy_switch = 2.5                           # s
        time_ideal, summary_ideal = delivery_timing.get_timing(tramp_e, x, y, tramp_w, energy_switch_time=False)
        total_ideal = time_ideal[-1]
        max_layers = np.floor((max_time - total_ideal)/default_energy_switch)
        print('Maximum theoretical energy layers in 2 minutes: {}'.format(max_layers))

        summary = summary.replace('Summary (s):', 'Original (s):')
        summary_adapted = summary_adapted.replace('Summary', 'Adapted')

        ax = fig.add_subplot(4, 2, 6)
        ax.plot(time, color='blue', alpha=0.5)
        ax.plot(time_adapted, color='red')
        accu_len = 0
        for layer_energy in unique_energies:
            time_lyr = time[tramp_e == layer_energy]
            time_adapted_lyr = time_adapted[tramp_e == layer_energy]
            add_rectangle(ax, accu_len, min(time_lyr), len(time_lyr), max(time_adapted_lyr)-min(time_lyr), 'green', 0.1)
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

        # FIGURE 1, PLOT 7 --------------------------------
        ax = fig.add_subplot(4, 2, 7)
        ax.hist(dw, 50, facecolor='green', normed=True, alpha=0.75)
        ax.set_xlabel("Weight change")
        ax.set_ylabel("Frequency")

        # FIGURE 1, PLOT 8 --------------------------------
        positions_keys = [str(i) for i in np.round(unique_energies, 1)]
        boxprops = dict(linewidth=0.6)
        whiskerprops = dict(linestyle='-', linewidth=0.5)
        medianprops = dict(linewidth=0.6)
        meanprops = dict(markersize=2)
        capprops = dict(linewidth=0.6)
        flierprops = dict(marker='o', markersize=1, markerfacecolor='black')
        box_kwargs = dict(return_type='dict', flierprops=flierprops, medianprops=medianprops,
                          boxprops=boxprops, whiskerprops=whiskerprops, patch_artist=True,
                          capprops=capprops, showmeans=True)

        # Create pandas df to plot boxplot
        weight_layers = list()
        for i,layer_energy in enumerate(unique_energies):
            weight_layers.append(dw[tramp_e == layer_energy])
        df = pd.DataFrame(weight_layers)
        df = df.transpose()
        df.columns = positions_keys

        ax = fig.add_subplot(4, 2, 8)
        bp = df.boxplot(ax=ax, meanprops=meanprops, **box_kwargs)
        ax.set_xlabel("Energy layer", fontsize=7)
        ax.set_ylabel(r'$\Delta w$', fontsize=7)
        ax.grid(False)
        ax.grid(color='k', linestyle=':', linewidth=0.5, alpha=0.25, axis='y', zorder=1)
        ax.set_xticklabels(labels=ax.get_xticklabels(), rotation=30, horizontalalignment='center')

        colors = ['#d9544d', '#3778bf', '#7bb274']
        [item.set_color('black') for item in bp['boxes']]
        [item.set_facecolor(colors[i%len(colors)]) for i,item in enumerate(bp['boxes'])]
        [item.set_color('black') for item in bp['whiskers']]
        [item.set_color('black') for item in bp['medians']]
        [item.set_markerfacecolor(colors[(i-1)%len(colors)]) for i,item in enumerate(bp['means'])]


        pp.savefig(fig, bbox_inches='tight')
        fig.clf()

        # FIGURE 2, PLOT 1 --------------------------------
        if spots_layer:
            max_ncols = 5
            ncols = min(number_layers, max_ncols)
            max_nrows = 5
            nrows = min(int(np.ceil(number_layers / ncols)), max_nrows)
            max_plots_page = nrows*ncols
            fig.suptitle('Spot movement per layer. Beam: {}. Layers: {}-{}'.format(
                               tramp_num, 0, min(max_plots_page-1, number_layers-1)))
            for layer_num, layer_energy in enumerate(unique_energies):
                subplot_num = layer_num % max_plots_page + 1
                ax_layers = fig.add_subplot(nrows, ncols, subplot_num)
                ax_layers.tick_params(labelsize=4)
                if subplot_num % ncols == 1:
                    ax_layers.set_ylabel('Y (mm)', fontsize=6)
                if int((subplot_num-1) / ncols) == nrows-1:
                    ax_layers.set_xlabel('X (mm)', fontsize=6)
                bool_mask = tramp_e == layer_energy
                temp_tramp_x = tramp_x[bool_mask]
                temp_tramp_y = tramp_y[bool_mask]
                temp_x = x[bool_mask]
                temp_y = y[bool_mask]
                temp_vx = temp_x - temp_tramp_x
                temp_vy = temp_y - temp_tramp_y
                temp_d = np.sqrt(temp_vx * temp_vx + temp_vy * temp_vy)

                for i in range(len(temp_tramp_x)):
                    if temp_d[i] > 0.001:
                        ax_layers.arrow(0.1 * temp_vx[i] + temp_tramp_x[i], 0.1 * temp_vy[i] + temp_tramp_y[i],
                                        0.8 * temp_vx[i], 0.8 * temp_vy[i],
                                        fc='k', ec='k', alpha=0.25, zorder=0)
                ax_layers.scatter(temp_tramp_x, temp_tramp_y, s=10, linewidth=0.5, zorder=1,
                                  edgecolors='black', facecolors='')
                ax_layers.scatter(temp_x, temp_y, s=10, linewidth=0.5, zorder=2, edgecolors='k',
                                  facecolors='blue')

                ax_layers.annotate(str(layer_energy) + ' MeV', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=5)
                if layer_num % max_plots_page + 1 == 0 or layer_num+1 == number_layers:
                    pp.savefig(fig, bbox_inches='tight')
                    fig.clf()
                    fig.suptitle('Spot movement per layer. Beam: {}. Layers: {}-{}'.format(
                        tramp_num, layer_num, min(layer_num + max_plots_page-1, number_layers-1)))


def main(argv):
    # Define the top-level parser
    parser = argparse.ArgumentParser(description='Analysis of adaptation')
    # Parse module:
    parser.add_argument('--vf',     help='File with vector field values and locations', required=True)
    parser.add_argument('--shifts', help='File with vector energy and locations shifts', required=True)
    parser.add_argument('--ctv',    help='MHA file containing the CTV structure')
    parser.add_argument('--tramps', nargs='+', help='MHA file containing the CTV structure', required=True)
    parser.add_argument('--outdir', help='Directory to output analysis', default='./')
    parser.add_argument('--spots_layer', help='Neglect layer-by-layer plotting of spot position shifts', default=False)
    parser.add_argument('--outfile', help='Report file name. It will be forced to be a pdf.',
                        default='adapt_report.pdf')

    args = parser.parse_args(argv)

    if not args.outfile.endswith('.pdf'):
        args.outfile += '.pdf'
    outfile = args.outdir + "/" + args.outfile
    pp = PdfPages(outfile)
    analize_vf(args.vf, pp)
    analize_tramp(args.shifts, args.tramps, args.spots_layer, pp)
    pp.close()


if __name__ == "__main__":
    main(sys.argv[1:])
