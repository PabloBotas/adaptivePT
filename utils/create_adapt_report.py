#!/usr/bin/python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from  more_itertools import unique_everseen
from mpl_toolkits.mplot3d.axes3d import Axes3D


# import matplotlib.gridspec as gridspec

# import collections

# import dvh_analysis
# import tester
# import reader
# import utils


def analize_vf(vf_file, outdir):
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

    # fig = plt.figure(figsize=(10, 8))

    ax = plt.subplot(2, 4, 1)
    Y,X = np.histogram(d, nbins)
    cm = plt.cm.get_cmap('YlOrRd')
    C = [cm(((i-X.min())/(X.max()-X.min()))) for i in X]
    ax.bar(X[:-1], Y, color=C, width=X[1]-X[0], alpha=0.75)

    ax.set_xlabel('Vector size (mm)', fontsize=11)
    ax.tick_params(labelsize=6)

    ax = plt.subplot(2, 4, 2, projection='polar')
    b, _, _ = ax.hist(ang_x, nangles, histtype='step', alpha = 1, color='r', fill=True, facecolor = 'r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle x', fontsize=12)
    ax.tick_params(labelsize=6)
    ax = plt.subplot(2, 4, 3, projection='polar')
    b, _, _ = ax.hist(ang_y, nangles, histtype='step', alpha = 1, color='r', fill=True, facecolor = 'r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle y', fontsize=12)
    ax.tick_params(labelsize=6)
    ax = plt.subplot(2, 4, 4, projection='polar')
    b, _, _ = ax.hist(ang_z, nangles, histtype='step', alpha = 1, color='r', fill=True, facecolor = 'r')
    ax.set_rticks(np.round(np.array([0.25, 0.5, 0.75, 1])*max(b)))
    ax.set_rmax(np.round(1.05*max(b)))
    ax.set_title('Angle z', fontsize=12)
    ax.tick_params(labelsize=6)
    
    ax = plt.subplot(2, 3, 4)
    ax.quiver(x, y, vx, vy, d, angles='xy', cmap=plt.cm.YlOrRd, pivot = 'tail', alpha = 0.75)
    ax.set_xlabel('pos x (mm)', fontsize=11)
    ax.set_ylabel('pos y (mm)', fontsize=11)
    ax.tick_params(labelsize=8)

    ax = plt.subplot(2, 3, 5)
    ax.quiver(x, z, vx, vz, d, angles='xy', cmap=plt.cm.YlOrRd, pivot = 'tail', alpha = 0.75)
    ax.set_xlabel('pos x (mm)', fontsize=11)
    ax.set_ylabel('pos z (mm)', fontsize=11)
    ax.tick_params(labelsize=8)
   
    ax = plt.subplot(2, 3, 6)
    ax.quiver(y, z, vy, vz, d, angles='xy', cmap=plt.cm.YlOrRd, pivot = 'tail', alpha = 0.75)
    ax.set_xlabel('pos y (mm)', fontsize=11)
    ax.set_ylabel('pos z (mm)', fontsize=11)
    ax.tick_params(labelsize=8)

    outfile = outdir + "/vector_field_analysis.pdf"
    plt.savefig(outfile, bbox_inches = 'tight')

def analize_tramp(shifts_file, tramp_files, outdir):
    r = np.genfromtxt(shifts_file, skip_header=1, delimiter = ' ',
                      names = ['e', 'x', 'y', 'z', 'beamid', 'spotid']).T
    all_de = r['e']/1e6
    all_vx = r['x']
    all_vy = r['y']
    beamid = r['beamid']
    spotid = r['spotid']

    nbins = 200
    nangles = 360
    cm = plt.cm.get_cmap('YlOrRd')

    for i, f in enumerate(tramp_files, start=0):
        r = np.genfromtxt(f, skip_header=12, names = ['e', 'x', 'y', 'w']).T
        e = r['e']
        x = r['x']
        y = r['y']
        w = r['w']

        de = all_de[beamid == i]
        vx = all_vx[beamid == i]
        vy = all_vy[beamid == i]
        d = np.sqrt(vx*vx + vy*vy)
        ang_x = np.arccos(vx/d)

        fig = plt.figure(figsize=(10, 8))

        ax = plt.subplot(2, 1, 1)
        accu_len = 0
        for j in np.unique(e)[::-1]:
            layer = (e+de)[e == j]
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
                    (accu_len, 0.95*j),  # (x,y)
                    len(layer),  # width
                    0.1*j,  # height
                    color='green',
                    alpha=0.075
                )
            )
            seg_x = [accu_len, accu_len+len(layer)]
            ax.plot(seg_x, [j, j], 'k', alpha=0.5)
            accu_len += len(layer)
        ax.plot(e+de, 'ro',  markersize=2)
        ax.set_xlabel('Spot number', fontsize=11)
        ax.set_ylabel('Energy (MeV)', fontsize=11)

        ax = plt.subplot(2, 1, 2)
        accu_len = 0
        for j in unique_everseen(y):
            layer = (y + vy)[y == j]
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
                    (accu_len, 0.95 * j),  # (x,y)
                    len(layer),  # width
                    0.1 * j,  # height
                    color='green',
                    alpha=0.075
                )
            )
            seg_x = [accu_len, accu_len + len(layer)]
            ax.plot(seg_x, [j, j], 'k', alpha=0.5)
            accu_len += len(layer)
        ax.plot(y + vy, color='blue', alpha=0.5)
        ax.plot(y + vy, 'ro', markersize=2)
        ax.set_xlabel('Spot number', fontsize=11)
        ax.set_ylabel('Slow dimension pos (mm)', fontsize=11)
        ax.annotate('Original energy layers = ' + str(len(np.unique(e))) + '\nAdapted energy layers = ' + str(len(np.unique(de))),
                    xy=(10, 10), xycoords='figure pixels')

        # ax = plt.subplot(1, 2, 2, projection='3d')
        # ax.view_init(10, 20)
        # phase_space_d = np.sqrt(vx*vx+vy*vy+de*de)
        # ax.quiver(x, y, e, vx, vy, de, phase_space_d, cmap=plt.cm.YlOrRd)


        filename = outdir + "/tramp_analysis_" + str(i) + ".pdf"
        plt.savefig(filename, bbox_inches = 'tight')
        plt.close()

def main(argv):
    ## Define the top-level parser
    parser = argparse.ArgumentParser(description='Analysis of adaptation')
    
    ## Parse split_tramps module:
    parser.add_argument('--vf',     help='File with vector field values and locations', required=True)
    parser.add_argument('--shifts', help='File with vector energy and locations shifts', required=True)
    parser.add_argument('--ctv',    help='MHA file containing the CTV structure', required=True)
    parser.add_argument('--tramps', nargs='+', help='MHA file containing the CTV structure', required=True)
    parser.add_argument('--outdir', help='Directory to output analysis', default='./')
    
    args = parser.parse_args()

    analize_vf(args.vf, args.outdir)
    analize_tramp(args.shifts, args.tramps, args.outdir)

if __name__=="__main__":
    main(sys.argv[1:])