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


def add_rectangle(ax, x, y, width, height, color, alpha):
    ax.add_patch(
        patches.Rectangle(
            (x, y), width, height,
            facecolor=color, alpha=alpha,
            edgecolor='black', linewidth=0.1
        )
    )


def analize_vf(vf_file, pp):
    r = np.genfromtxt(vf_file, skip_header=1, delimiter=' ',
                      names=['vx', 'vy', 'vz', 'x', 'y', 'z', 'beamid', 'spotid']).T
    vx = np.array(r['vx'])
    vy = np.array(r['vy'])
    vz = np.array(r['vz'])
    x = np.array(r['x'])
    y = np.array(r['y'])
    z = np.array(r['z'])
    # beamid = r['beamid']
    # spotid = r['spotid']

    npoints = len(x)

    d = np.sqrt(vx*vx + vy*vy + vz*vz)
    ang_x = np.array([np.arccos(vx[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])
    ang_y = np.array([np.arccos(vy[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])
    ang_z = np.array([np.arccos(vz[i]/d[i]) if d[i] != 0 else 0 for i in range(npoints)])

    nbins = 25
    nangles = 360

    fig = plt.figure()
    fig.suptitle('Vector field analysis at endpoints')

    ax = fig.add_subplot(2, 4, 1)
    bins_y, bins_x = np.histogram(d, nbins)
    cm = plt.cm.get_cmap('rainbow')
    hist_colors = [cm(((i-bins_x.min())/(bins_x.max()-bins_x.min()))) for i in bins_x]
    ax.bar(bins_x[:-1], bins_y, color=hist_colors, width=bins_x[1]-bins_x[0], alpha=0.75)
    ax.set_xlabel('Vector size (mm)', fontsize=8)

    ax = fig.add_subplot(2, 4, 2, projection='polar')
    b, _, _ = ax.hist(ang_x, nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle x', fontsize=11)
    ax = fig.add_subplot(2, 4, 3, projection='polar')
    b, _, _ = ax.hist(ang_y, nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle y', fontsize=11)
    ax = fig.add_subplot(2, 4, 4, projection='polar')
    b, _, _ = ax.hist(ang_z, nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle z', fontsize=11)

    ax = fig.add_subplot(2, 3, 4)
    dummy = vx if vx.any() or vy.any() else np.full((npoints, 1), 0.00000001)
    ax.quiver(x, y, dummy, vy, d, angles='xy', scale_units='xy', scale=1,
              cmap=plt.cm.get_cmap('rainbow'), pivot='tail', alpha=0.75)
    ax.set_xlabel('pos x (mm)', fontsize=8)
    ax.set_ylabel('pos y (mm)', fontsize=8)

    ax = fig.add_subplot(2, 3, 5)
    dummy = vy if vy.any() or vz.any() else np.full((npoints, 1), 0.00000001)
    ax.quiver(y, z, dummy, vz, d, angles='xy', scale_units='xy', scale=1,
              cmap=plt.cm.get_cmap('rainbow'), pivot='tail', alpha=0.75)
    ax.set_xlabel('pos y (mm)', fontsize=8)
    ax.set_ylabel('pos z (mm)', fontsize=8)

    ax = fig.add_subplot(2, 3, 6)
    dummy = vz if vz.any() or vx.any() else np.full((npoints, 1), 0.00000001)
    ax.quiver(z, x, dummy, vx, d, angles='xy', scale_units='xy', scale=1,
              cmap=plt.cm.get_cmap('rainbow'), pivot='tail', alpha=0.75)
    ax.set_xlabel('pos z (mm)', fontsize=8)
    ax.set_ylabel('pos x (mm)', fontsize=8)

    pp.savefig(fig, bbox_inches='tight')


def analize_tramp(shifts_file, tramp_files, spots_layer, pp):
    r = np.genfromtxt(shifts_file, skip_header=1, delimiter=' ',
                      names=['e', 'x', 'y', 'z', 'd', 'beamid', 'spotid']).T
    all_de = np.array(r['e']/1e6)
    all_x = 10*np.array(r['x'])
    all_y = 10*np.array(r['y'])
    # all_z  = np.array(r['z'])
    # all_d  = np.array(r['d'])
    beamid = np.array(r['beamid'])
    # spotid = r['spotid']

    for tramp_num, tramp_file in enumerate(tramp_files, start=0):
        tramp_r = np.genfromtxt(tramp_file, skip_header=12, names=['e', 'x', 'y', 'w']).T
        tramp_e = np.array(tramp_r['e'])
        tramp_x = np.array(tramp_r['x'])
        tramp_y = np.array(tramp_r['y'])
        tramp_w = np.array(tramp_r['w'])

        de = np.array(all_de[beamid == tramp_num])
        x = np.array(all_x[beamid == tramp_num])
        y = np.array(all_y[beamid == tramp_num])
        vx = x - tramp_x
        vy = y - tramp_y
        d = np.sqrt(vx*vx + vy*vy)

        unique_energies = np.unique(tramp_e)[::-1]
        number_layers = len(unique_energies)
        layer_id = np.array([int(np.squeeze(np.where(unique_energies == i))) for i in tramp_e])

        fig = plt.figure()
        fig.suptitle('Tramp adaptations analysis. Beam: {}'.format(tramp_num))

        # FIGURE 1, PLOT 1 --------------------------------
        ax = fig.add_subplot(3, 2, 1)
        accu_len = 0
        for i, layer_energy in enumerate(unique_energies):
            bool_mask = tramp_e == layer_energy
            layer = (tramp_e+de)[bool_mask]
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
        ax.plot(np.sort(tramp_e+de)[::-1], linestyle='None', marker='o', color='black',
                markersize=1, alpha=0.4, markeredgecolor='black', markeredgewidth=0.1)
        ax.plot(tramp_e+de, linestyle='None', color='red', marker='o', markersize=1.75, markeredgewidth=0.25)
        ax.set_xlabel('Spot number', fontsize=7)
        ax.set_ylabel('Energy (MeV)', fontsize=7)
        ax.annotate('Number spots = ' + str(len(tramp_e)) + '\nOriginal energy layers  = ' +
                    str(number_layers) + '\nAdapted energy layers = ' + str(len(np.unique(tramp_e+de))),
                    xy=(5, min(min(tramp_e+de), min(tramp_e))), fontsize=6,
                    textcoords='axes fraction', xytext=(0.04, 0.04))

        # FIGURE 1, PLOT 2 --------------------------------
        ax = fig.add_subplot(3, 2, 2)
        nbins = 25
        ax.hist(de, nbins)
        ax.set_xlabel('Shift size (MeV)', fontsize=7)

        # FIGURE 1, PLOT 3 --------------------------------
        ax = fig.add_subplot(3, 2, 3)
        for i in range(len(tramp_x)):
            if d[i] > 0.001:
                ax.arrow(0.1 * vx[i] + tramp_x[i], 0.1 * vy[i] + tramp_y[i],
                         0.8 * vx[i], 0.8 * vy[i],
                         fc='k', ec='k', alpha=0.25, zorder=0)
        cmap = plt.cm.get_cmap('hsv')
        scat_colors = np.array(cmap(layer_id / layer_id.max()))
        ax.scatter(tramp_x, tramp_y, s=10, linewidth=0.5, zorder=1,
                   edgecolors=scat_colors - np.array([0, 0, 0, 0.75]), facecolors='')
        ax.scatter(x, y, s=10, linewidth=0.5, zorder=2, edgecolors='k',
                   facecolors=scat_colors - np.array([0, 0, 0, 0]))
        ax.set_xlabel('X (mm)', fontsize=7)
        ax.set_ylabel('Y (mm)', fontsize=7)

        # FIGURE 1, PLOT 4 --------------------------------
        ax = fig.add_subplot(3, 2, 4)
        nbins = 25
        if d.max()-d.min() < 0.1:
            ax.hist(d, bins=nbins, range=(d.min-0.05, d.max+0.05))
        else:
            ax.hist(d, bins=nbins)
        ax.set_xlabel('Shift size (mm)', fontsize=7)

        # FIGURE 1, PLOT 5 --------------------------------
        ax = fig.add_subplot(3, 2, 5)
        ax.step(range(len(tramp_y)), tramp_y, color='black', alpha=0.5)
        ax.step(range(len(tramp_y)), y, color='blue', alpha=0.5)
        ax.plot(y, linestyle='None', linewidth=0.5, color='red', marker='o', markersize=2)
        ax.set_xlabel('Spot number', fontsize=7)
        ax.set_ylabel('Slow dimension pos (mm)', fontsize=7)
        ax.annotate('Number spots = ' + str(len(tramp_e)) + '\n' +
                    'Original slow layers = ' + str(len(np.unique(tramp_y))) + '\n' +
                    'Adapted slow layers = ' + str(len(np.unique(y))),
                    xy=(5, min(min(y), min(tramp_y))), fontsize=6,
                    textcoords='axes fraction', xytext=(0.04, 0.04))

        time, summary = delivery_timing.get_timing(tramp_e, tramp_x, tramp_y, tramp_w)
        time_adapted, summary_adapted = delivery_timing.get_timing(tramp_e+de, x, y, tramp_w)
        max_time = 120. if 120. > time[-1] else 1.5*time[-1]  # s
        default_energy_switch = 2.5                           # s
        time_ideal, summary_ideal = delivery_timing.get_timing(tramp_e, x, y, tramp_w, energy_switch_time=False)
        total_ideal = time_ideal[-1]
        max_layers = np.floor((max_time - total_ideal)/default_energy_switch)
        print('Maximum theoretical energy layers: {}'.format(max_layers))

        summary = summary.replace('Summary (s):', 'Original (s):')
        summary_adapted = summary_adapted.replace('Summary', 'Adapted')

        # FIGURE 1, PLOT 6 --------------------------------
        ax = fig.add_subplot(3, 2, 6)
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
