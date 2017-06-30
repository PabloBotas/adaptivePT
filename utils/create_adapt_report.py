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
    vx = r['vx']
    vy = r['vy']
    vz = r['vz']
    x = r['x']
    y = r['y']
    z = r['z']
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
    all_de = r['e']/1e6
    all_vx = r['x']
    all_vy = r['y']
    beamid = r['beamid']
    # spotid = r['spotid']

    # nbins = 200
    # nangles = 360
    # cm = plt.cm.get_cmap('YlOrRd')

    figures = list()

    for tramp_num, tramp_file in enumerate(tramp_files, start=0):
        r = np.genfromtxt(tramp_file, skip_header=12, names=['e', 'x', 'y', 'w']).T
        e = r['e']
        x = r['x']
        y = r['y']
        w = r['w']

        de = all_de[beamid == tramp_num]
        vx = all_vx[beamid == tramp_num]
        vy = all_vy[beamid == tramp_num]
        # d = np.sqrt(vx*vx + vy*vy)
        # ang_x = np.arccos(vx/d)

        fig = plt.figure()
        fig.suptitle('Tramp adaptations analysis. Beam: {}'.format(tramp_num))

        ax = fig.add_subplot(2, 2, 1)
        accu_len = 0
        for layer_energy in np.unique(e)[::-1]:
            layer = (e+de)[e == layer_energy]
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
        ax.plot(e+de, 'ro',  markersize=2)
        ax.set_xlabel('Spot number', fontsize=8)
        ax.set_ylabel('Energy (MeV)', fontsize=8)
        ax.annotate('Number spots = ' + str(len(e)) + '\nOriginal energy layers  = ' +
                    str(len(np.unique(e))) + '\nAdapted energy layers = ' + str(len(np.unique(e+de))),
                    xy=(5, min(min(e+de), min(e))), fontsize=6,
                    textcoords='axes fraction', xytext=(0.04, 0.04))
        ax.tick_params(labelsize=8)

        ax = fig.add_subplot(2, 2, 3)
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

        ax = fig.add_subplot(1, 2, 2)
        ax.plot(time, color='blue', alpha=0.5)
        ax.plot(time_adapted, color='red')
        accu_len = 0
        for layer_energy in np.unique(e)[::-1]:
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
                    textcoords='axes fraction', xytext=(0.04, 0.81),
                    xy=(5, min(time)), fontsize=6, bbox=bbox_props)
        ax.annotate(summary_adapted, family='monospace',
                    textcoords='axes fraction', xytext=(0.04, 0.615),
                    xy=(5, min(time)), fontsize=6, bbox=bbox_props_adapt)

        # filename = outdir + "/tramp_analysis_" + str(tramp_num) + ".pdf"
        # plt.savefig(filename, bbox_inches='tight')
        # plt.close()

        figures.append(fig)
    return figures


def main(argv):
    # Define the top-level parser
    parser = argparse.ArgumentParser(description='Analysis of adaptation')
    
    # Parse split_tramps module:
    parser.add_argument('--vf',     help='File with vector field values and locations', required=True)
    parser.add_argument('--shifts', help='File with vector energy and locations shifts', required=True)
    parser.add_argument('--ctv',    help='MHA file containing the CTV structure', required=True)
    parser.add_argument('--tramps', nargs='+', help='MHA file containing the CTV structure', required=True)
    parser.add_argument('--outdir', help='Directory to output analysis', default='./')
    parser.add_argument('--split',  help='If output PDF should be splitted', default=False)

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
