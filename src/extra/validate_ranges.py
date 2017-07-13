#!/usr/bin/python3

import sys
import argparse
import numpy as np
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy import ndimage

mpl.rcParams['axes.labelpad'] = 2
mpl.rcParams['axes.formatter.useoffset'] = False
mpl.rcParams['font.size'] = 6
mpl.rcParams['figure.titlesize'] = 14
mpl.rcParams['figure.dpi'] = 250


def read_info(topas_files, tramp_files):
    info = list()

    # Read Topas files --------
    couch_string = 'd:Rt/beam/PatientSupportAngle='
    gantry_string = 'd:Rt/beam/Gantry='
    dose_nx_string = 'i:Rt/dose/Rows='
    dose_ny_string = 'i:Rt/dose/Columns='
    dose_nz_string = 'i:Rt/dose/NumberOfSlicesZ='
    dose_dx_string = 'd:Rt/dose/PixelSpacing0='
    dose_dy_string = 'd:Rt/dose/PixelSpacing1='
    dose_dz_string = 'd:Rt/dose/DoseGridSliceThickness='
    ct_nx_string = 'i:Rt/CT/Rows='
    ct_ny_string = 'i:Rt/CT/Columns='
    ct_nz_string = 'i:Rt/CT/NumberOfSlices='
    ct_dx_string = 'd:Rt/CT/PixelSpacing0='
    ct_dy_string = 'd:Rt/CT/PixelSpacing1='
    ct_dz_string = 'dv:Rt/CT/SliceThicknessSpacing='

    itera = topas_files if isinstance(topas_files, list) else [topas_files]
    for i, file in enumerate(itera):
        info.append(dict())
        for line in open(file, 'r'):
            if couch_string in line:
                info[i]['couch'] = float(line.split(' ')[1])
            if gantry_string in line:
                info[i]['gantry'] = float(line.split(' ')[1])
            if dose_nx_string in line:
                info[i]['dose_nx'] = int(line.split(' ')[1])
            if dose_ny_string in line:
                info[i]['dose_ny'] = int(line.split(' ')[1])
            if dose_nz_string in line:
                info[i]['dose_nz'] = int(line.split(' ')[1])
            if dose_dx_string in line:
                info[i]['dose_dx'] = float(line.split(' ')[1])
            if dose_dy_string in line:
                info[i]['dose_dy'] = float(line.split(' ')[1])
            if dose_dz_string in line:
                info[i]['dose_dz'] = float(line.split(' ')[1])
            if ct_nx_string in line:
                info[i]['ct_nx'] = int(line.split(' ')[1])
            if ct_ny_string in line:
                info[i]['ct_ny'] = int(line.split(' ')[1])
            if ct_nz_string in line:
                info[i]['ct_nz'] = int(line.split(' ')[1])
            if ct_dx_string in line:
                info[i]['ct_dx'] = float(line.split(' ')[1])
            if ct_dy_string in line:
                info[i]['ct_dy'] = float(line.split(' ')[1])
            if ct_dz_string in line:
                info[i]['ct_dz'] = float(line.split(' ')[2])
        info[i]['dose_n'] = info[i]['dose_nx'] * info[i]['dose_ny'] * info[i]['dose_nz']
        info[i]['ct_n'] = info[i]['ct_nx'] * info[i]['ct_ny'] * info[i]['ct_nz']

    # Read Tramp files --------
    itera = tramp_files if isinstance(tramp_files, list) else [tramp_files]
    for i, file in enumerate(itera):
        e = list()
        x = list()
        y = list()
        w = list()
        for line in open(file, 'r'):
            if '#' in line:
                continue
            e_, x_, y_, w_ = line.split('\t')
            e.append(float(e_))
            x.append(float(x_))
            y.append(float(y_))
            w.append(float(w_))
        info[i]['e'] = np.array(e)
        info[i]['x'] = np.array(x)
        info[i]['y'] = np.array(y)
        info[i]['w'] = np.array(w)

    return info


def set_info_dims(info, l):
    for i, d in enumerate(info):
        if d['dose_n'] == l:
            d['nx'] = d['dose_nx']
            d['ny'] = d['dose_ny']
            d['nz'] = d['dose_nz']
            d['dx'] = d['dose_dx']
            d['dy'] = d['dose_dy']
            d['dz'] = d['dose_dz']
            if i == 0:
                print('The files contain dose grid dimensions...')
        else:
            d['nx'] = d['ct_nx']
            d['ny'] = d['ct_ny']
            d['nz'] = d['ct_nz']
            d['dx'] = d['ct_dx']
            d['dy'] = d['ct_dy']
            d['dz'] = d['ct_dz']
            if i == 0:
                print('The files contain CT dimensions...')
    return info


def get_interpolated_index(y, array, shift=0):
    i = shift + np.argmin(np.abs(array[shift:] - y))
    x0 = i - 1 if array[i] > y else i
    y0 = array[x0]
    y1 = array[x0 + 1]
    # print((y-y0)/(y1-y0), i, shift+(y-y0)/(y1-y0), shift)
    return x0 + (y-y0)/(y1-y0)
    # return x0


def get_ranges(mc_list, rays_list, info, nspots, plot_profiles, outfile):
    itera_mc = mc_list if isinstance(mc_list, list) else [mc_list]
    itera_ray = rays_list if isinstance(rays_list, list) else [rays_list]
    assert len(itera_mc) == len(itera_ray)
    spots_per_beam = int(len(itera_mc)/len(info))

    nbeams = len(info)
    rows_per_figure = 4 if 4 < nspots else nspots
    ranges = dict()
    ranges['ray'] = list()
    ranges['max'] = list()
    ranges['r90'] = list()
    ranges['r80'] = list()
    ranges['r70'] = list()
    ranges['r50'] = list()
    ranges['factor'] = list()
    if plot_profiles:
        pp = PdfPages("{}.pdf".format(outfile))
        fig = plt.figure()

    for ibeam in range(nbeams):
        print("Working on beam {}:".format(ibeam))
        if plot_profiles:
            fig.suptitle("Beam {}".format(ibeam))

        for idose in range(nspots):
            if plot_profiles and idose % rows_per_figure == 0 and idose > 0:
                pp.savefig(fig, bbox_inches='tight')
                fig.clf()
            mc_data = np.fromfile(itera_mc[ibeam*spots_per_beam+idose], dtype=np.float32)
            ray_data = np.fromfile(itera_ray[ibeam*spots_per_beam+idose], dtype=np.float32)
            assert len(mc_data) == len(ray_data)

            if idose == 0 and ibeam == 0:
                info = set_info_dims(info, len(mc_data))
                ranges['x'] = np.arange(info[ibeam]['nx']) * info[ibeam]['dx']

            if nspots < 10:
                print("    Spot {} ...".format(idose))
            elif idose % 10 == 0:
                print("    Spot {} ...".format(idose))
            mc_data = mc_data.reshape((info[ibeam]['nz'], info[ibeam]['ny'], info[ibeam]['nx']))
            ray_data = ray_data.reshape((info[ibeam]['nz'], info[ibeam]['ny'], info[ibeam]['nx']))
            mc_data = mc_data.sum(axis=0)
            ray_data = ray_data.any(axis=0)

            y, x = np.where(ray_data > 0.)
            x = x[::-1]
            y = y[::-1]
            pol = np.poly1d(np.polyfit(x, y, 1))
            angle = np.arctan(pol[1])

            temp = mc_data
            mc_data = ndimage.rotate(mc_data, 180./np.pi*angle, reshape=False, order=0)
            ray_data = ndimage.rotate(ray_data, 180./np.pi*angle, reshape=False, order=0)
            if plot_profiles:
                ax = fig.add_subplot(rows_per_figure, 3, 3*(idose % rows_per_figure) + 1)
                ax.imshow(temp, aspect='auto', cmap='YlOrRd_r', origin='lower', alpha=0.5)
                ax.imshow(ray_data, aspect='auto', interpolation='none', cmap='binary', origin='lower', alpha=0.5)
                ax.imshow(mc_data, aspect='auto', cmap='YlOrRd_r', origin='lower', alpha=0.5)
                ax.annotate("Spot {}".format(idose), xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                            xycoords=ax.yaxis.label, textcoords='offset points',
                            size='large', ha='right', va='center', rotation=90)

            mc_dd = mc_data.sum(axis=0)
            ray_dd = ray_data.any(axis=0)
            if plot_profiles:
                ax = fig.add_subplot(rows_per_figure, 3, 3*(idose % rows_per_figure) + 2)
                ax.plot(mc_dd, color='black')
                ax2 = ax.twinx()
                ax2.plot(ray_dd, linewidth=0.2, color='blue')
                ax2.yaxis.set_ticklabels([])

            mc_prof = mc_data.sum(axis=1)
            ray_prof = ray_data.any(axis=1)
            if plot_profiles:
                ax = fig.add_subplot(rows_per_figure, 3, 3*(idose % rows_per_figure) + 3)
                ax.plot(mc_prof, color='black')
                ax2 = ax.twinx()
                ax2.plot(ray_prof, linewidth=0.2, color='blue')
                ax2.yaxis.set_ticklabels([])

            # Find values
            m = np.max(mc_dd)
            ranges['factor'].append(info[ibeam]['dx']/np.cos(angle))
            ind = np.max(np.where(mc_dd == m))
            ranges['ray'].append(np.max(np.nonzero(ray_dd)))
            ranges['max'].append(ind)
            ranges['r90'].append(get_interpolated_index(0.9*m, mc_dd, shift=ind))
            ranges['r80'].append(get_interpolated_index(0.8*m, mc_dd, shift=ind))
            ranges['r70'].append(get_interpolated_index(0.7*m, mc_dd, shift=ind))
            ranges['r50'].append(get_interpolated_index(0.5*m, mc_dd, shift=ind))

        if plot_profiles:
            pp.savefig(fig, bbox_inches='tight')
            fig.clf()
    if plot_profiles:
        pp.close()

    ranges['factor'] = np.array(ranges['factor'])
    ranges['max'] = np.array(ranges['max'])
    ranges['r90'] = np.array(ranges['r90'])
    ranges['r80'] = np.array(ranges['r80'])
    ranges['r70'] = np.array(ranges['r70'])
    ranges['r50'] = np.array(ranges['r50'])

    energies = list()
    for i in info:
        energies.append(i['e'][:nspots])
    energies = np.array(energies).flatten()
    dev_max = np.array(100 * (ranges['ray'] - ranges['max']) / ranges['max'])
    dev_r90 = np.array(100 * (ranges['ray'] - ranges['r90']) / ranges['r90'])
    dev_r80 = np.array(100 * (ranges['ray'] - ranges['r80']) / ranges['r80'])
    dev_r70 = np.array(100 * (ranges['ray'] - ranges['r70']) / ranges['r70'])
    dev_r50 = np.array(100 * (ranges['ray'] - ranges['r50']) / ranges['r50'])

    with open('{}.txt'.format(outfile), 'w', newline='') as file:
        for i in range(len(energies)):
            file.write('{} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(
                energies[i], ranges['factor'][i], ranges['ray'][i], ranges['max'][i],
                ranges['r90'][i], ranges['r80'][i], ranges['r70'][i], ranges['r50'][i],
                dev_max[i], dev_r90[i], dev_r80[i], dev_r70[i], dev_r50[i]
            ))


def analyze_ranges(infile, outfile):

    ranges = dict()
    ranges['factor'] = list()
    ranges['ray'] = list()
    ranges['max'] = list()
    ranges['r90'] = list()
    ranges['r80'] = list()
    ranges['r70'] = list()
    ranges['r50'] = list()
    energies = list()
    dev_max = list()
    dev_r90 = list()
    dev_r80 = list()
    dev_r70 = list()
    dev_r50 = list()
    with open(infile, 'r') as file:
        for line in file:
            e, f, rr, rm, r90, r80, r70, r50, dm, d90, d80, d70, d50 = line.split()
            energies.append(float(e))
            ranges['factor'].append(float(f))
            ranges['ray'].append(float(rr))
            ranges['max'].append(float(rm))
            ranges['r90'].append(float(r90))
            ranges['r80'].append(float(r80))
            ranges['r70'].append(float(r70))
            ranges['r50'].append(float(r50))
            dev_max.append(float(dm))
            dev_r90.append(float(d90))
            dev_r80.append(float(d80))
            dev_r70.append(float(d70))
            dev_r50.append(float(d50))

    ranges['factor'] = np.array(ranges['factor'])
    ranges['ray'] = np.array(ranges['ray'])
    ranges['max'] = np.array(ranges['max'])
    ranges['r90'] = np.array(ranges['r90'])
    ranges['r80'] = np.array(ranges['r80'])
    ranges['r70'] = np.array(ranges['r70'])
    ranges['r50'] = np.array(ranges['r50'])

    fig = plt.figure()
    fig.suptitle("Range analysis")
    ax = fig.add_subplot(1, 1, 1)
    # ax.plot(energies, dev_max, marker='o', color='black', linestyle='None', label='Max')
    ax.plot(energies, dev_r90, marker='o', color='blue',  linestyle='None', label='R90')
    ax.plot(energies, dev_r80, marker='o', color='red',   linestyle='None', label='R80')
    ax.plot(energies, dev_r70, marker='o', color='green', linestyle='None', label='R70')
    # ax.plot(energies, dev_r50, marker='o', color='gray',  linestyle='None', label='R50')
    ax.legend()

    j = 0
    for i in np.arange(len(dev_r80)):
        if abs(dev_r80[i]) > 2:
            ax.annotate('{:0.2f} : Beam {} - Spot {}'.format(
                            dev_r80[i], int(i/(len(dev_r80)/2)), int(i % (len(dev_r80)/2))
                        ),
                        xy=(0, 0), xytext=(1.08, 1-j*0.05),
                        xycoords='axes fraction', textcoords='axes fraction',
                        size='small', ha='center', va='center')
            j += 1

    # pol_max = np.poly1d(np.polyfit(energies, dev_max, 1))
    pol_r90 = np.poly1d(np.polyfit(energies, dev_r90, 1))
    pol_r80 = np.poly1d(np.polyfit(energies, dev_r80, 1))
    pol_r70 = np.poly1d(np.polyfit(energies, dev_r70, 1))
    # pol_r50 = np.poly1d(np.polyfit(energies, dev_r50, 1))
    # ax.plot(energies, pol_max(energies), color='black')
    ax.plot(energies, pol_r90(energies), color='blue')
    ax.plot(energies, pol_r80(energies), color='red')
    ax.plot(energies, pol_r70(energies), color='green')
    # ax.plot(energies, pol_r50(energies), color='gray')

    pp = PdfPages("{}.pdf".format(outfile))
    pp.savefig(fig, bbox_inches='tight')
    pp.close()


def main(argv):
    # Define the top-level parser
    parser = argparse.ArgumentParser(description='Analysis of adaptation')
    # Parse module:
    parser.add_argument('--mc', nargs='+', help='Binary file with ray tracing of each spot')
    parser.add_argument('--rays', nargs='+', help='Binary file with dose of each spot')
    parser.add_argument('--tramps', nargs='+', help='Tramp files')
    parser.add_argument('--nspots', help='Number of spots to analyze per tramp file', default=20, type=int)
    parser.add_argument('--info', nargs='+', help='Topas files with the angles and dimensions')
    parser.add_argument('--input', help='Input data for analysis')
    parser.add_argument('--outdir', help='Output directory', default='./')
    parser.add_argument('--outfiles', help='Output files basename', required=True)

    args = parser.parse_args(argv)

    if args.input is None and (not args.mc or not args.rays or not args.info or not args.tramps):
        print("ERROR! Wrong arguments. Please, check help")
        parser.print_help()
        exit(1)

    args.outdir = args.outdir[:-1] if args.outdir.endswith('/') else args.outdir
    basefile = "{}/{}".format(args.outdir, args.outfiles)
    basefile = basefile if basefile.rfind('.') < basefile.rfind('/') else basefile[:basefile.rfind('.')]
    if args.input is None:
        info = read_info(args.info, args.tramps)
        get_ranges(args.mc, args.rays, info, args.nspots, True, basefile)
    else:
        analyze_ranges(args.input, basefile)

if __name__ == "__main__":
    main(sys.argv[1:])
