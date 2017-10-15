#!/usr/bin/python3

import sys
import argparse
import numpy as np
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy import stats
import pandas as pd

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

            ## Fit to get the angle and rotate the beam for analysis
            y, x = np.where(ray_data > 0.)
            x = x[::-1]
            y = y[::-1]
            pol = np.poly1d(np.polyfit(x, y, 1))
            angle = np.arctan(pol[1])

            temp = mc_data
            mc_data = ndimage.rotate(mc_data, 180./np.pi*angle, reshape=False, order=0)
            ray_data = ndimage.rotate(ray_data, 180./np.pi*angle, reshape=False, order=0)
            y, x = np.where(ray_data > 0.)
            pol = np.poly1d(np.polyfit(x, y, 1))
            if plot_profiles:
                ax = fig.add_subplot(rows_per_figure, 3, 3*(idose % rows_per_figure) + 1)
                ax.imshow(temp, aspect='auto', cmap='YlOrRd_r', origin='lower', alpha=0.5)
                ax.plot(x, pol(x), color='black')
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


def my_linfit(x, y):
    z = np.polyfit(x,y,1)
    p = np.poly1d(z)
    p_y = z[0]*x + z[1]
    y_err = y - p_y
    # p_x = np.arange(np.min(x),np.max(x)+1,1)
    p_x = np.unique(x)
    # now calculate confidence intervals for new test x-series
    mean_x = np.mean(x) # mean of x
    n = len(x)          # number of samples in original fit
    t = stats.t.ppf(1-0.025, len(x))  # appropriate t value (two tailed 95%)
    s_err = np.sum(np.power(y_err,2))  # sum of the squares of the residuals
    confs = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
                ((np.sum(np.power(x,2)))-n*(np.power(mean_x,2))))))
    # now predict y based on test x-values
    p_y = z[0]*p_x+z[1]
    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)
    return p, upper, lower


def analyze_ranges(infile, outfile):
    plt.style.use('ggplot')
    colors = ['#d9544d', '#3778bf', '#7bb274']
    plt.locator_params(nbins=10)
    df = pd.read_csv(infile, sep=" ", header=None)
    df.columns = ['energy', 'factor', 'ray', 'max',
                  'r90', 'r80', 'r70', 'r50', 'dev_max',
                  'dev_r90', 'dev_r80', 'dev_r70', 'dev_r50']
    df['energy_round'] = df.apply(lambda row: int(10*row.energy+0.5)/10, axis=1)
    
    fig, ax = plt.subplots(nrows=2, ncols=2)
    ## PLOT boxplots
    boxprops = dict(linewidth=0.6)
    whiskerprops = dict(linestyle='-', linewidth=0.5)
    medianprops = dict(linewidth=0.6)
    meanprops_1 = dict(markersize=3)
    meanprops_2 = dict(markersize=5)
    capprops = dict(linewidth=0.6)
    flierprops = dict(marker='o', markersize=1, markerfacecolor='black')
    box_kwargs = dict(return_type='dict', flierprops=flierprops, medianprops=medianprops,
                      boxprops=boxprops, whiskerprops=whiskerprops, patch_artist=True,
                      capprops=capprops, showmeans=True)
    positions = np.unique(df['energy'])
    bp0 = df.boxplot(['dev_r90'], by=['energy_round'], ax=ax[0,0], widths=1.5, positions=positions, meanprops=meanprops_1, **box_kwargs)   
    bp1 = df.boxplot(['dev_r80'], by=['energy_round'], ax=ax[0,1], widths=1.5, positions=positions, meanprops=meanprops_1, **box_kwargs)
    bp2 = df.boxplot(['dev_r70'], by=['energy_round'], ax=ax[1,0], widths=1.5, positions=positions, meanprops=meanprops_1, **box_kwargs)
    bp3 = df.boxplot(['dev_r90','dev_r80','dev_r70'], ax=ax[1,1], meanprops=meanprops_2, **box_kwargs)

    ## Set colors (f***ing pandas boxprops do not work)
    [item.set_color('black') for item in bp3['boxes']]
    [item.set_facecolor(colors[i]) for i,item in enumerate(bp3['boxes'])]
    [item.set_color('black') for item in bp3['whiskers']]
    [item.set_color('black') for item in bp3['medians']]
    [item.set_markerfacecolor(colors[i-1]) for i,item in enumerate(bp3['means'])]
    for i, bp in enumerate([bp0, bp1, bp2]):
        [item.set_color('black') for item in bp[list(bp.keys())[0]]['boxes']]
        [item.set_facecolor(colors[i]) for item in bp[list(bp.keys())[0]]['boxes']]
        [item.set_color('black') for item in bp[list(bp.keys())[0]]['whiskers']]
        [item.set_color('black') for item in bp[list(bp.keys())[0]]['medians']]
        [item.set_markerfacecolor(colors[i-1]) for item in bp[list(bp.keys())[0]]['means']]

    ## PLOT TREND LINES AND CI
    for i,(var,axid) in enumerate(zip(['dev_r90', 'dev_r80', 'dev_r70'], [(0,0),(0,1),(1,0)])):
        fit, upper, lower = my_linfit(df['energy'].values, df[var].values)
        ax[axid].plot(df['energy'].values, fit(df['energy'].values), color=colors[i], alpha=0.9)
        ax[axid].fill_between(np.unique(df['energy'].values), lower, upper, color=colors[i], alpha=0.4)
        ax[axid].annotate('{}'.format(fit),
                          xy=(0, 0), xytext=(0.5, 0.95),
                          xycoords='axes fraction', textcoords='axes fraction',
                          size='small', ha='center', va='center')
    ## PLOT AVERAGES AND SIGMAS
    anno_kwargs = dict(xy=(0, 0),
                       xycoords='axes fraction', textcoords='axes fraction',
                       size='small', ha='center', va='center')
    ax[(1,1)].annotate(r'${:.2f} \pm {:.2f}$'.format(df['dev_r90'].mean(), df['dev_r90'].std()), xytext=(0.165, 0.95), **anno_kwargs)
    ax[(1,1)].annotate(r'${:.2f} \pm {:.2f}$'.format(df['dev_r80'].mean(), df['dev_r80'].std()), xytext=(0.500, 0.95), **anno_kwargs)
    ax[(1,1)].annotate(r'${:.2f} \pm {:.2f}$'.format(df['dev_r70'].mean(), df['dev_r70'].std()), xytext=(0.825, 0.95), **anno_kwargs)

    fig.suptitle("")
    fig.tight_layout()
    
    titles = ['R90', 'R80', 'R70', 'Comparison']
    for i,pair in enumerate([(0,0),(0,1),(1,0),(1,1)]):
        ax[pair].xaxis.label.set_fontsize(10)
        ax[pair].yaxis.label.set_fontsize(10)
        ax[pair].set(title=titles[i], xlabel='Energy (MeV)', ylabel='Deviation (%)')
        ax[pair].title.set_fontsize(10)
        rotation = 0 if pair == (1,1) else 90
        ax[pair].set_xticklabels(ax[pair].get_xticklabels(), rotation=rotation, size=7)
        ax[pair].set_yticklabels(ax[pair].get_yticklabels(), size=7)

    ## FIND AND ADD OUTLIERS INFORMATION
    j = 0
    for i in df.index:
        if abs(df['dev_r80'][i]) > 2:
            ax[0,1].annotate('{:0.2f} : Beam {} - Spot {}'.format(
                    df['dev_r80'][i], int(i/(len(df.index)/2)), int(i % (len(df.index)/2))
                ),
                xy=(0, 0), xytext=(1.24, 1-j*0.05),
                xycoords='axes fraction', textcoords='axes fraction',
                size='small', ha='center', va='center')
            j += 1

    ## SAVE TO FILE
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
    parser.add_argument('--input', help='Processed input data for analysis')
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
