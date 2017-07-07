#!/usr/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from more_itertools import unique_everseen
import delivery_timing
from matplotlib.backends.backend_pdf import PdfPages

def analize_vf(vf_file):
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

    d = np.sqrt(vx*vx + vy*vy + vz*vz)
    ang_x = np.arccos(vx/d)
    ang_y = np.arccos(vy/d)
    ang_z = np.arccos(vz/d)

    nbins = 25
    nangles = 360

    fig = plt.figure()
    fig.suptitle('Vector field analysis at endpoints')

    ax = fig.add_subplot(2, 4, 1)
    bins_y, bins_x = np.histogram(d, nbins)
    cm = plt.cm.get_cmap('YlOrRd')
    colors = [cm(((i-bins_x.min())/(bins_x.max()-bins_x.min()))) for i in bins_x]
    ax.bar(bins_x[:-1], bins_y, color=colors, width=bins_x[1]-bins_x[0], alpha=0.75)

    ax.set_xlabel('Vector size (mm)', fontsize=11)
    ax.tick_params(labelsize=6)

    ax = fig.add_subplot(2, 4, 2, projection='polar')
    b, _, _ = ax.hist(ang_x, nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle x', fontsize=12)
    ax.tick_params(labelsize=6)
    ax = fig.add_subplot(2, 4, 3, projection='polar')
    b, _, _ = ax.hist(ang_y, nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle y', fontsize=12)
    ax.tick_params(labelsize=6)
    ax = fig.add_subplot(2, 4, 4, projection='polar')
    b, _, _ = ax.hist(ang_z, nangles, histtype='step', alpha=1, color='r', fill=True, facecolor='r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle z', fontsize=12)
    ax.tick_params(labelsize=6)
    
    ax = fig.add_subplot(2, 3, 4)
    ax.quiver(x, y, vx, vy, d, angles='xy', cmap=plt.cm.YlOrRd, pivot='tail', alpha=0.75)
    ax.set_xlabel('pos x (mm)', fontsize=11)
    ax.set_ylabel('pos y (mm)', fontsize=11)
    ax.tick_params(labelsize=8)

    ax = fig.add_subplot(2, 3, 5)
    ax.quiver(x, z, vx, vz, d, angles='xy', cmap=plt.cm.YlOrRd, pivot='tail', alpha=0.75)
    ax.set_xlabel('pos x (mm)', fontsize=11)
    ax.set_ylabel('pos z (mm)', fontsize=11)
    ax.tick_params(labelsize=8)
   
    ax = fig.add_subplot(2, 3, 6)
    ax.quiver(y, z, vy, vz, d, angles='xy', cmap=plt.cm.YlOrRd, pivot='tail', alpha=0.75)
    ax.set_xlabel('pos y (mm)', fontsize=11)
    ax.set_ylabel('pos z (mm)', fontsize=11)
    ax.tick_params(labelsize=8)

    # outfile = outdir + "/vector_field_analysis.pdf"
    # plt.savefig(outfile, bbox_inches='tight')
    return fig


def analize_tramp(shifts_file, tramp_files):
    r = np.genfromtxt(shifts_file, skip_header=1, delimiter=' ',
                      names=['e', 'x', 'y', 'z', 'beamid', 'spotid']).T
    all_de = np.array(r['e']/1e6)
    all_vx = np.array(r['x'])
    all_vy = np.array(r['y'])
    beamid = np.array(r['beamid'])
    # spotid = r['spotid']

    # nbins = 200
    # nangles = 360
    # cm = plt.cm.get_cmap('YlOrRd')

    figures = list()

    for tramp_num, tramp_file in enumerate(tramp_files, start=0):
        r = np.genfromtxt(tramp_file, skip_header=12, names=['e', 'x', 'y', 'w']).T
        e = np.array(r['e'])
        x = np.array(r['x'])
        y = np.array(r['y'])
        w = np.array(r['w'])

        de = np.array(all_de[beamid == tramp_num])
        vx = np.array(all_vx[beamid == tramp_num])
        vy = np.array(all_vy[beamid == tramp_num])
        # d = np.sqrt(vx*vx + vy*vy)
        # ang_x = np.arccos(vx/d)

        unique_energies = np.unique(e)[::-1]
        number_layers = len(unique_energies)

        fig_general = plt.figure()
        fig_general.suptitle('Tramp adaptations analysis. Beam: {}'.format(tramp_num))

        ax = fig_general.add_subplot(2, 2, 1)
        accu_len = 0
        for layer_energy in unique_energies:
            bool_mask = e == layer_energy
            layer = (e+de)[bool_mask]
            average = np.mean(layer)
            pos_text_y = 1.01*max(layer) if max(layer) > 1.05*layer_energy else 1.01*1.05*layer_energy
            ax.annotate("{:.2f} %".format(
                100*(average-layer_energy)/layer_energy), xy=(accu_len, pos_text_y), fontsize=4)
            ax.add_patch(
                patches.Rectangle(
                    (accu_len, min(layer)),  # (x,y)
                    len(layer),  # width
                    max(layer)-min(layer),  # height
                    color='blue',
                    alpha=0.1
                )
            )
            ax.add_patch(
                patches.Rectangle(
                    (accu_len, 0.95*layer_energy),  # (x,y)
                    len(layer),  # width
                    0.1*layer_energy,  # height
                    color='green',
                    alpha=0.075
                )
            )
            ax.plot([accu_len, accu_len+len(layer)], [layer_energy, layer_energy], 'k', alpha=0.5)
            accu_len += len(layer)
        ax.plot(np.sort(e+de)[::-1], 'ko', markersize=1, alpha=0.2)
        ax.plot(e+de, 'ro', markersize=2)
        ax.set_xlabel('Spot number', fontsize=8)
        ax.set_ylabel('Energy (MeV)', fontsize=8)
        ax.annotate('Number spots = ' + str(len(e)) + '\nOriginal energy layers  = ' +
                    str(number_layers) + '\nAdapted energy layers = ' + str(len(np.unique(e+de))),
                    xy=(5, min(min(e+de), min(e))), fontsize=6,
                    textcoords='axes fraction', xytext=(0.04, 0.04))
        ax.tick_params(labelsize=8)

        ax = fig_general.add_subplot(2, 2, 3)
        accu_len = 0
        previous_layer_energy = 0
        for layer_energy in unique_everseen(y):
            layer = (y + vy)[y == layer_energy]
            ax.add_patch(
                patches.Rectangle(
                    (accu_len, min(layer)),  # (x,y)
                    len(layer),  # width
                    max(layer) - min(layer),  # height
                    color='blue',
                    alpha=0.1
                )
            )
            ax.add_patch(
                patches.Rectangle(
                    (accu_len, 0.95 * layer_energy),  # (x,y)
                    len(layer),  # width
                    0.1 * layer_energy,  # height
                    color='green',
                    alpha=0.075
                )
            )
            ax.plot([accu_len, accu_len + len(layer)], [layer_energy, layer_energy], 'k', alpha=0.5)
            ax.plot([accu_len, accu_len], [previous_layer_energy, layer_energy],
                    linestyle=':', color='black', alpha=0.25)
            previous_layer_energy = layer_energy
            accu_len += len(layer)
        ax.plot(y + vy, color='blue', alpha=0.5)
        ax.plot(y + vy, 'ro', markersize=2)
        ax.set_xlabel('Spot number', fontsize=8)
        ax.set_ylabel('Slow dimension pos (mm)', fontsize=8)
        ax.annotate('Number spots = ' + str(len(e)) + '\nOriginal slow direction layers  = ' +
                    str(len(np.unique(y))) + '\nAdapted slow direction layers = ' + str(len(np.unique(y + vy))),
                    xy=(5, min(min(y + vy), min(y))), fontsize=6,
                    textcoords='axes fraction', xytext=(0.04, 0.04))
        ax.tick_params(labelsize=8)

        time, summary = delivery_timing.get_timing(e, x, y, w)
        time_adapted, summary_adapted = delivery_timing.get_timing(e+de, x+vx, y+vy, w)
        summary = summary.replace('Summary (s):', 'Original (s):')
        summary_adapted = summary_adapted.replace('Summary', 'Adapted')

        ax = fig_general.add_subplot(2, 2, 2)
        ax.plot(time, color='blue', alpha=0.5)
        ax.plot(time_adapted, color='red')
        accu_len = 0
        for layer_energy in unique_energies:
            time_lyr = time[e == layer_energy]
            time_adapted_lyr = time_adapted[e == layer_energy]
            ax.add_patch(
                patches.Rectangle(
                    (accu_len, min(time_lyr)),  # (x,y)
                    len(time_lyr),  # width
                    max(time_adapted_lyr) - min(time_lyr),  # height
                    color='green',
                    alpha=0.1
                )
            )
            accu_len += len(time_lyr)
        ax.set_xlabel('Spot number', fontsize=8)
        ax.set_ylabel('Time (s)', fontsize=8)
        ax.tick_params(labelsize=8)
        bbox_props = dict(boxstyle="Round", facecolor="blue", alpha=0.1)
        bbox_props_adapt = dict(boxstyle="Round", facecolor="red", alpha=0.1)
        ax.annotate(summary, family='monospace',
                    textcoords='axes fraction', xytext=(0.03, 0.77),
                    xy=(5, min(time)), fontsize=4, bbox=bbox_props)
        ax.annotate(summary_adapted, family='monospace',
                    textcoords='axes fraction', xytext=(0.32, 0.77),
                    xy=(5, min(time)), fontsize=4, bbox=bbox_props_adapt)

        ax = fig_general.add_subplot(2, 2, 4)
        ax.scatter(x, y, facecolors='none', edgecolors='black', alpha=0.5)
        for i in range(len(x)):
            ax.arrow(0.1 * vx[i] + x[i], 0.1 * vy[i] + y[i],
                     0.80 * vx[i], 0.80 * vy[i],
                     fc='k', ec='k', alpha=0.5)
        ax.scatter(x + vx, y + vy, edgecolor='')
        ax.tick_params(labelsize=8)
        ax.set_ylabel('Y (mm)', fontsize=8)
        ax.set_xlabel('X (mm)', fontsize=8)

        figures.append(fig_general)

        # filename = outdir + "/tramp_analysis_" + str(tramp_num) + ".pdf"
        # plt.savefig(filename, bbox_inches='tight')
        # plt.close()

        # PLOT SPOTS SHIFTS PER LAYER
        max_ncols = 5
        ncols = min(number_layers, max_ncols)
        max_nrows = 5
        nrows = min(int(np.ceil(number_layers / ncols)), max_nrows)
        max_plots_page = nrows*ncols
        fig_spots = plt.figure()
        fig_spots.suptitle('Spot movement per layer. Beam: {}. Layers: {}-{}'.format(
                           tramp_num, 0, min(max_plots_page-1, number_layers-1)))
        for layer_num, layer_energy in enumerate(unique_energies):
            subplot_num = layer_num % max_plots_page + 1
            ax_layers = fig_spots.add_subplot(nrows, ncols, subplot_num)
            ax_layers.tick_params(labelsize=4)
            if subplot_num % ncols == 1:
                ax_layers.set_ylabel('Y (mm)', fontsize=6)
            if int((subplot_num-1) / ncols) == nrows-1:
                ax_layers.set_xlabel('X (mm)', fontsize=6)
            bool_mask = e == layer_energy
            pos_x = x[bool_mask]
            pos_y = y[bool_mask]
            shift_x = vx[bool_mask]
            shift_y = vy[bool_mask]
            ax_layers.scatter(pos_x, pos_y, facecolors='none', edgecolors='black', alpha=0.5)
            for i in range(len(pos_x)):
                ax_layers.arrow(0.1*shift_x[i]+pos_x[i], 0.1*shift_y[i]+pos_y[i],
                                0.80*shift_x[i], 0.80*shift_y[i],
                                fc='k', ec='k', alpha=0.5)
            ax_layers.scatter(pos_x+shift_x, pos_y+shift_y, edgecolor='')
            ax_layers.annotate(str(layer_energy) + ' MeV', xy=(0.05, 0.05), xycoords='axes fraction', fontsize=5)
            if layer_num % max_plots_page + 1 == 0 or layer_num+1 == number_layers:
                figures.append(fig_spots)
                fig_spots = plt.figure()
                fig_spots.suptitle('Spot movement per layer. Beam: {}. Layers: {}-{}'.format(
                    tramp_num, layer_num, min(layer_num + max_plots_page-1, number_layers-1)))

    return figures


def main(argv):
    # Define the top-level parser
    parser = argparse.ArgumentParser(description='Analysis of adaptation')
    
    # Parse split_tramps module:
    parser.add_argument('--vf',     help='File with vector field values and locations', required=True)
    parser.add_argument('--shifts', help='File with vector energy and locations shifts', required=True)
    parser.add_argument('--ctv',    help='MHA file containing the CTV structure', required=False)
    parser.add_argument('--tramps', nargs='+', help='MHA file containing the CTV structure', required=True)
    parser.add_argument('--outdir', help='Directory to output analysis', default='./')
    parser.add_argument('--split',  help='If output PDF should be split', default=False)

    args = parser.parse_args()

    fig_vf = analize_vf(args.vf)
    figs_tramp = analize_tramp(args.shifts, args.tramps)

    if args.split:
        file = args.outdir + "/vector_field_analysis.pdf"
        plt.savefig(file, bbox_inches='tight')
        for i, fig in enumerate(figs_tramp):
            file = args.outdir + "/tramp_analysis_" + str(i) + ".pdf"
            plt.savefig(file, bbox_inches='tight')
    else:
        pp = PdfPages(args.outdir + "/adapt_report.pdf")
        pp.savefig(fig_vf, bbox_inches='tight')
        for i in figs_tramp:
            pp.savefig(i, bbox_inches='tight')
        pp.close()

if __name__ == "__main__":
    main(sys.argv[1:])
