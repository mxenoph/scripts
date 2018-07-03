#!/usr/bin/env python

# Import modules# {{{
# Needed for division to return float
from __future__ import division
import os, sys, argparse, re
import warnings
from time import gmtime, strftime
import fnmatch
import gffutils
import numpy as np
import pybedtools
import metaseq
from pylab import *
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
from collections import OrderedDict
warnings.simplefilter("error")
# }}}

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()
# Required arguments
parser.add_argument('-p', '--path', metavar = "path", type = str, required = True, help = 'Path to search for npz files')
# Optional arguments
parser.add_argument('-m', '--matching', type=str, required = True, help = 'Pattern by which to retrieve files from the path. e.g. "*3KO*"',
        # default pattern set to all so I don't have to write an if else condition for defining files
        default = "*")
parser.add_argument('-s', '--subsets', type=str, required = False, help = 'Tab delimited file with gene_id, group columns')

args = parser.parse_args()

if not args.path.endswith(os.sep):
    args.path = args.path + os.sep

plot_path = args.path + 'plots/'

# Create plot directory if it doesn't exist
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

#array_order = ['5000-gene_start.npz', 'genes-filtered.npz', 'gene_end-5000.npz']
array_order = ['2000-gene_start.npz', 'genes-filtered.npz', 'gene_end-500.npz']
tmp = [ os.path.join(args.path, f) for f in os.listdir(args.path) if f.endswith(".npz") and fnmatch.fnmatch(f, args.matching) ]
#tmp = [ f for f in tmp if f.endswith('.5000-gene_start.npz') or f.endswith('.genes-filtered.npz') or f.endswith('.gene_end-5000.npz') ]
tmp = [ f for f in tmp if f.endswith('.2000-gene_start.npz') or f.endswith('.genes-filtered.npz') or f.endswith('.gene_end-500.npz') ]
assert len(tmp) == 3, \
        'Input arrays provided are not 3.'

# Not the best way but only way I can ensure the order of the arrays in the tuple is gene start, gene, gene end
files = []
for m in array_order:
    for f in tmp:
        if f.endswith(m):
            files = files + [f,]
print files
#}}}

def get_N_HexCol(N=15):# {{{
    import colorsys
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#' + "".join(map(lambda x: chr(x).encode('hex'),rgb)))
    return hex_out
# }}}

def format_default_axes(figure, sample, feature_name = "Genes", ticks = [-2000, 0, 500], tick_labels = ['-2Kb', 'TSS', '500bp'], alignment = ['right', 'center', 'right']):# {{{
    figure.line_axes.axis('tight')
    for ax in [figure.line_axes, figure.array_axes]:
        ax.axvline(0, color='k', linestyle=':')
        if 'TTS' in tick_labels:
            ax.axvline(ticks[tick_labels.index('TTS')], color='k', linestyle=':')
        
    # Label axes
    figure.array_axes.set_title(sample)
    figure.array_axes.set_ylabel("Genes")

    figure.line_axes.set_ylabel("Average enrichment (IP - input)\n RPMMR")
    figure.line_axes.set_xlabel(feature_name)

    figure.cax.set_ylabel("RPMMR")
    
    # Remove ticks from heatmap x axis
    for tick in figure.array_axes.get_xticklabels():
        tick.set_visible(False)

    # Set ticks and labels on the average profile plot
    figure.line_axes.set_xticks(ticks)
    figure.line_axes.set_xticklabels(tick_labels)
    for txt, ha in zip(figure.line_axes.get_xticklabels(), alignment):
        txt.set_horizontalalignment(ha)

    figure.subplots_adjust(left = 0.2, right = 0.8, hspace = 0.05, top = 0.95)

    # Positioning the colobar axes took some tweaking:
    figure.cax.set_position([0.75, .15, .05, .10])
    figure.set_facecolor('w')

    return(figure)
# }}}

def format_subsets_axes(figure, by, order, ncols):# {{{
    # Add average plot per subsets
    print "In format_subsets_axes subsets"
    metaseq.plotutils.add_labels_to_subsets(
            figure.array_axes,
            subset_by = by,
            subset_order = order)
    figure.line_axes.legend(loc ='best', frameon = False)

    box = figure.line_axes.get_position()
    figure.line_axes.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    figure.line_axes.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.3), ncol = ncols, frameon = False)

    return(figure)
# }}}

def plot_subsets_separate(array, arguments, pp, sample, feature_name, ticks = [-2000, 0, 500], tick_labels = ['-2Kb', 'TSS', '500bp'], alignment = ['right', 'center', 'right']):# {{{
    inds = []
    for x in arguments['subset_order']:
        subset_ind = np.nonzero(arguments['subset_by'] == x)[0]
        subset_sort_by = arguments['sort_by'][subset_ind]
        subset_argsort_by = np.argsort(subset_sort_by)
        inds.append(subset_ind[subset_argsort_by])
    ind = np.concatenate(inds)

    from matplotlib import pyplot as plt
    plt.rcParams['font.size'] = 11
    plt.rcParams['legend.scatterpoints'] = 1
    plt.rcParams['legend.fontsize'] = 10
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for subset_ind, label, _lkw, _fkw in zip(inds, arguments['subset_order'], arguments['line_kwargs'], arguments['fill_kwargs']):
        metaseq.plotutils.ci_plot(arguments['x'],
                array[subset_ind],
                ax = ax,
                line_kwargs = _lkw,
                fill_kwargs = _fkw)

    # http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    ax.set_title(sample)

    # Add a vertical line at TSS/gene-start
    ax.axvline(0, color='k', linestyle=':')
    if 'TTS' in tick_labels:
        ax.axvline(ticks[tick_labels.index('TTS')], color='k', linestyle=':')

    ax.set_xlabel(feature_name)
    ax.set_ylabel('Average enrichment (IP-input) RPMMR')

    ax.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.1), ncol = len(arguments['fill_kwargs']), frameon = False)
    
    # Set ticks and labels on the average profile plot
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    for txt, ha in zip(ax.get_xticklabels(), alignment):
        txt.set_horizontalalignment(ha)

    pp.savefig(fig)
    plt.close(fig)
# }}}

def format_expression_axes(figure, expression, quantiles = None, e_type = "log(FPKM)", ticks = [-2000, 0, 500], tick_labels = ['-2Kb', 'TSS', '500bp'], alignment = ['right', 'center', 'right']):# {{{
    print "Formatting expression panel"
    # Axes limit tweaks
    figure.strip_axes.axis('tight')
    figure.strip_axes.axis(xmin = 0)

    # Remove all edges but the bottom x-axis.
    figure.strip_axes.yaxis.set_visible(False)
    figure.strip_axes.spines['top'].set_visible(False)
    figure.strip_axes.spines['right'].set_visible(False)
    figure.strip_axes.spines['left'].set_visible(False)
    figure.strip_axes.xaxis.set_ticks_position('bottom')

    figure.strip_axes.set_xlabel(e_type)
    figure.strip_axes.xaxis.set_major_locator(
        matplotlib.ticker.MaxNLocator(nbins = 5, prune = None))
    
    # Add a vertical line at TSS/gene-start
    for ax in [figure.line_axes, figure.array_axes]:
        ax.axvline(0, color='k', linestyle=':')
        if 'TTS' in tick_labels:
            ax.axvline(ticks[tick_labels.index('TTS')], color='k', linestyle=':')

    # Set ticks and labels on the average profile plot
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels)
    for txt, ha in zip(ax.get_xticklabels(), alignment):
        txt.set_horizontalalignment(ha)

    # Plot expression data in the right-most panel
    figure.strip_axes.fill_betweenx(
            x1 = sorted(expression),
            y = np.arange(len(expression)),
            x2 = 0,
            color = '.5')

    # If information about silent, Q1-Q4 genes is given add a subplot to the figure# {{{
    if quantiles is not None:
        qlabel = {
                    0: 'silent',
                    1: '0-25%', # Coming from rscript 1 is Q1, 2 is Q2 etc
                    2: '25-50%',
                    3: '50-75%',
                    4: '75-100%'}
        # Create linearly spaced colors along a colormap, but with some padding at
        # the beginning to avoid too-light colors.
        quantile_colors = matplotlib.cm.YlOrBr(np.linspace(.5, 1, 4 + 1))
        uquantiles = sorted(list(set(quantiles)))
        
        # Create another axes along the bottom.
        quantile_ax = plt.subplot(figure.gs[2, 0])
        # Plot the quantile lines as different colors
        for q, color in zip(uquantiles, quantile_colors)[::-1]:
            ind = (quantiles == q).values
            metaseq.plotutils.ci_plot(
                    np.arange(700),
                    narray[ind, :],
                    ax = quantile_ax,
                    line_kwargs = dict(color = color, label = qlabel[q], alpha = 0.5),
                    fill_kwargs = dict(color = color, label = qlabel[q], alpha = 0.3)
                )
    
        # Axes configure on the bottom axes.
#        bottom_ax.set_xticks(ticks)
#        bottom_ax.set_xticklabels(ticklabels)
#        for txt, ha in zip(bottom_ax.get_xticklabels(), alignment):
#            txt.set_horizontalalignment(ha)
#        bottom_ax.axis('tight')
#        bottom_ax.axvline(100, **border_style)
#        bottom_ax.axvline(600, **border_style)
#        bottom_ax.legend(loc='best', fontsize=8, frameon=False, labelspacing=0.4)
#        bottom_ax.set_ylabel('Enrichment\n(RPMMR)')# }}}
    
    return(figure)
# }}}

# match_array_features: fset is a dataframe with indeces the gene_ids# {{{
def match_array_features(array, features, fset):
    fset_order = fset.index
    # Filtering based on original order as reindexing in line 173 will result in all features in array been represented in fset
    keep = [i for i, gene in enumerate(features) if gene.name in fset_order]
    # subset array and features
    s_array = array[keep]
    s_features = features.filter(lambda gene: gene.name in fset.index).saveas()

    fset_set = fset.reindex_to(s_features, attribute = 'gene_id')
    print "Features originally in array: %s \nFeatures in file provided: %s.\nFeatures removed from array: %s" % (array.shape[0], len(fset), array.shape[0] - len(fset_set))
    # If I do not write the feature_ids to variable and retrieve them on the fly in the tuple, it takes a long time
    feature_ids = [g.name for g in features]
    tmp = [x for x in fset.index if x not in feature_ids]
    print "Features removed from subsetting file %s:" % len(tmp)

    return {'s_array':s_array, 's_features':s_features, 's_fset':fset_set, 'fset_order': fset_order}
# }}}

# Plotting metagene# {{{
def plot_metagene(narray, features, label = None, fsets = None):

    # Default params# {{{
    arguments = {
            'vmin':5, 'vmax':95, 'percentile':True,
            'figsize':(5, 8),
            'line_kwargs': dict(color = 'k'),
            'fill_kwargs': dict(color='k', alpha=0.4),
            # For metagenes filter based on max value on each row of the metagene array
            # because genes with somewhat signal at the TSS might not be ordered after high signal
            # TSSs if the gene body signal drives the average down
            'sort_by': narray.max(axis=1)
            }

    plt.rcParams['font.family'] = 'Sans'
    plt.rcParams['font.size'] = 10

    # bins are 10bp long so to get the gene to look nice and long as many bins as we have I need to multiply the number
    # of bins with 10
    window = np.linspace(-5000, 0, 500)
    window = np.append(window, np.linspace(0,10000,1000))
    window = np.append(window, np.linspace(10000,15000,500))
    arguments['x'] = window

#    ticks, tick_labels, alignment = zip(*[
#        (-5000, '-5kb', 'right'),
#        (0, 'TSS', 'center'),
#        (2500, '', 'center'),
#        (5000, '50%', 'center'),
#        (7500, '', 'center'),
#        (10000, 'TTS', 'center'),
#        (15000, '5kb', 'left'),
#        ])
#    
    ticks, tick_labels, alignment = zip(*[
        (-2000, '-2kb', 'right'),
        (0, 'TSS', 'center'),
        (2500, '', 'center'),
        (5000, '50%', 'center'),
        (7500, '', 'center'),
        (10000, 'TTS', 'center'),
        (10500, '500b', 'left'),
        ])

    # }}}

    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages(plot_path + label + '-metagenes-profiles.pdf')
    
    print "Plotting without subsetting or reordering"
    fig = metaseq.plotutils.imshow(narray, features = features, **arguments)
    fig = format_default_axes(fig, sample = label + ' (' + str(narray.shape[0]) + ' features)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)
    pp.savefig(fig)
    plt.close(fig)

    if fsets is not None:
        assert type(fsets) is not 'metaseq.results_table.ResultsTable', \
                'fsets passed is not a table.'

        assert len([i for i in ['groups', 'expression', 'order'] if i in list(fsets.columns.values)]) != 0, \
                'Table does not contain any column that I know how to handle'

        # Making sure that fsets and array have the same dimensions# {{{
        if len(features) != len(fsets):
            print "From now on working with sub-setted array and features."
            s_array_and_features = match_array_features(array = narray, features = features, fset = fsets)
            
            narray = s_array_and_features['s_array']
            features = s_array_and_features['s_features']
            fsets = s_array_and_features['s_fset']
            fset_order = s_array_and_features['fset_order']

            # There is no point in plotting the heatmap again if the array is not sub-setted.
            # Reindexing is done in match_array_features(), but if dimensions are the same it's done in else condition 
            # before continuing with any grouping etc
            fig = metaseq.plotutils.imshow(narray, features = features, **arguments)
            fig = format_default_axes(fig, sample = label + ' (' + str(narray.shape[0]) + ' features)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)
            pp.savefig(fig)
            plt.close(fig)
        else:
            print "Re-indexing fsets provided."
            fsets = fsets.reindex_to(features, attribute = 'gene_id')

        # }}}

        if 'order' in list(fsets.columns.values):# {{{
            print "Order column provided...Subsetting and/or ordering"
            # http://nbviewer.jupyter.org/github/daler/metaseq/blob/master/doc/source/example_session.ipynb
            # matplotlib.imshow with the argument `origin="lower"`, which means the row in the plot at y=0
            # corresponds to the last row in the array (index=-1).
            arguments['sort_by'] = fsets.order
            o_fig = metaseq.plotutils.imshow(narray, **arguments)
            o_fig = format_default_axes(o_fig, sample = label + ' (' + str(narray.shape[0]) + ' features - ordered)', feature_name = "Metagenes")
            pp.savefig(o_fig)
            plt.close(o_fig)
            # }}}

        if 'groups' in list(fsets.columns.values): # {{{
            print "plot_metagene(): subsetting based on groups provided"
            groups = fsets.groups.unique()
            cls = np.zeros(len(narray)).astype('str')
            
            # The numpy arrays for each group label will have the same length as len(features) == len(fsets),
            # and same order because fsets was reindexed before.
            subset = []
            for g in groups:
                subset.append((str(g), (fsets.groups == g).values))

            for lab, ind in subset:
                cls[ind] = str(lab)
            # if group names are numbers then if 0 is a group asset will fail even if all features are assigned to a group
            assert sum(cls == '0.0') == 0

            s_arguments  = arguments
            # Saving groups and fsets to arguments for plotutils
            s_arguments['subset_by'] = cls
            s_arguments['subset_order'] = map(str, sorted(groups))

            s_arguments['line_kwargs'] = []
            s_arguments['fill_kwargs'] = []
            colors = get_N_HexCol(len(groups))
            
            for i, g in enumerate(sorted(groups)):
                s_arguments['line_kwargs'].append(dict(color=colors[i], label = g))
                s_arguments['fill_kwargs'].append(dict(color=colors[i], alpha = 0.3))
            
            s_fig = metaseq.plotutils.imshow(narray, **s_arguments)
            s_fig = format_default_axes(s_fig, sample = label + ' (' + str(narray.shape[0]) + ' features)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)
            s_fig = format_subsets_axes(s_fig, by = s_arguments['subset_by'], order = s_arguments['subset_order'], ncols = len(arguments['fill_kwargs']))

            pp.savefig(s_fig)
            plt.close(s_fig)
            
            plot_subsets_separate(array= narray, arguments = s_arguments, pp = pp, sample = label + ' (' + str(narray.shape[0]) + ' features)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)

            if 'order' in list(fsets.columns.values):
            # ordering within the cluster is done from bottom to top of heatmap e.g
            # e.g if ordering is the df on the left, on the heatmap from top to bottom will be
            # a 1 classA    classB c
            # b 2 classA    classB d
            # c 3 classB    classA a
            # d 4 classB    classA b
                arguments['sort_by'] = fsets.order
                
                s_fig = metaseq.plotutils.imshow(narray, **s_arguments)
                s_fig = format_default_axes(s_fig, sample = label + ' (' + str(narray.shape[0]) + ' features - ordered)', feature_name = "Metagenes")
                s_fig = format_subsets_axes(s_fig, by = s_arguments['subset_by'], order = s_arguments['subset_order'], ncols = len(s_arguments['fill_kwargs']))
                pp.savefig(s_fig)
                plt.close(s_fig)
        # }}}

        if 'expression' in list(fsets.columns.values) and False:# {{{
            print "Printing expression panel."
            e_arguments = arguments
            e_arguments['sort_by'] = narray.max(axis=1)

            # Params for the strip plot alongside the array
            e_arguments['subplot_params'] = dict(wspace=0.1, hspace=0.1, bottom=0.05)
            e_arguments['strip'] = True
            e_arguments['width_ratios'] = (4, 1)
            e_arguments['height_ratios'] = (4, 1, 1)
            
            e_fig = metaseq.plotutils.imshow(narray, **e_arguments)
            e_fig = format_default_axes(e_fig, sample = label + ' (' + str(narray.shape[0]) + ' features - ordered by expression)', feature_name = "Metagenes")
            e_fig = format_expression_axes(e_fig, expression = fsets.expression, ticks = ticks, tick_labels = tick_labels, alignment = alignment)
            pp.savefig(e_fig)
            plt.close(e_fig)
        # }}}

    pp.close()
    return {'narray': narray, 'features': features, 'fsets': fsets, 'fset_order':fset_order}
# }}}

def plot_array(array, features, pp, label = None, fsets = None):# {{{
    arguments = {
            'vmin':5, 'vmax':95, 'percentile':True,
            'figsize':(5, 8),
            'line_kwargs': dict(color = 'k'),
            'fill_kwargs': dict(color='k', alpha=0.4),
            'fill_kwargs': dict(color='k', alpha=0.4),
            'sort_by': array.mean(axis=1)
            }

    plt.rcParams['font.family'] = 'Sans'
    plt.rcParams['font.size'] = 10

    # Calculate number of bins depending on feature# {{{
    pattern = re.compile('(\d+)-(\w+)-(\d+)')
    # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
    matching = pattern.match(label)

    if matching:
        upstream, feature_type, downstream = matching.groups()
        # Set the bins such that counting is done for every 10bp window
        bins = (int(upstream) + int(downstream))/10
        window = np.linspace(-int(upstream), int(downstream), bins)
    elif re.compile('(\d+)-(\w+)').match(label) or re.compile('(\w+)-(\d+)').match(label):
        if re.compile('(\d+)-(\w+)').match(label):
            upstream, feature_type = re.compile('(\d+)-(\w+)').match(label).groups()
            bins = int(upstream)/10
            window = np.linspace(-int(upstream), 0, bins)
        else:
            feature_type, downstream = re.compile('(\w+)-(\d+)').match(label).groups()
            bins = int(downstream)/10
            window = np.linspace(0, int(downstream), bins)
    else:
        # For genes it only makes sense to count every 1% of gene
        bins = 1000
        window = np.linspace(0, 100, bins)
    #}}}

    arguments['x'] = window

    print "For checking how array %s should look before merging." % label
    fig = metaseq.plotutils.imshow(array, features = features, **arguments)
    fig = format_default_axes(fig, sample = label, feature_name = label, ticks = ticks, tick_labels = tick_labels, alignment = alignment)
    pp.savefig(fig)
    plt.close(fig)
    
    if fsets is not None:
        assert type(fsets) is not 'metaseq.results_table.ResultsTable', \
                'fsets passed is not a table.'

        assert len([i for i in ['groups', 'expression', 'order'] if i in list(fsets.columns.values)]) != 0, \
                'Table does not contain any column that I know how to handle'

        # Making sure that fsets and array have the same dimensions# {{{
        if len(features) != len(fsets):
            print "From now on working with sub-setted array and features."
            s_array_and_features = match_array_features(array = array, features = features, fset = fsets)
            
            array = s_array_and_features['s_array']
            features = s_array_and_features['s_features']
            fsets = s_array_and_features['s_fset']
            fset_order = s_array_and_features['fset_order']

            # There is no point in plotting the heatmap again if the array is not sub-setted.
            # Reindexing is done in match_array_features(), but if dimensions are the same it's done in else condition 
            # before continuing with any grouping etc
            fig = metaseq.plotutils.imshow(array, features = features, **arguments)
            fig = format_default_axes(fig, sample = label + ' (' + str(array.shape[0]) + ' features)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)
            pp.savefig(fig)
            plt.close(fig)
        else:
            print "Re-indexing fsets provided."
            fsets = fsets.reindex_to(features, attribute = 'gene_id')

        # }}}

        if 'order' in list(fsets.columns.values) and False:# {{{
            print "Order column provided...Subsetting and/or ordering"
            arguments['sort_by'] = fsets.order
            o_fig = metaseq.plotutils.imshow(array, **arguments)
            o_fig = format_default_axes(o_fig, sample = label + ' (' + str(array.shape[0]) + ' features - ordered)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)
            pp.savefig(o_fig)
            plt.close(o_fig)
            # }}}

        if 'groups' in list(fsets.columns.values): # {{{
            print "plot_metagene(): subsetting based on groups provided"
            groups = fsets.groups.unique()
            cls = np.zeros(len(array)).astype('str')
            
            # The numpy arrays for each group label will have the same length as len(features) == len(fsets),
            # and same order because fsets was reindexed before.
            subset = []
            for g in groups:
                subset.append((str(g), (fsets.groups == g).values))

            for lab, ind in subset:
                cls[ind] = str(lab)
            # if group names are numbers then if 0 is a group asset will fail even if all features are assigned to a group
            assert sum(cls == '0.0') == 0

            s_arguments  = arguments
            # Saving groups and fsets to arguments for plotutils
            s_arguments['subset_by'] = cls
            s_arguments['subset_order'] = map(str, sorted(groups))

            s_arguments['line_kwargs'] = []
            s_arguments['fill_kwargs'] = []
            colors = get_N_HexCol(len(groups))
            
            for i, g in enumerate(groups):
                s_arguments['line_kwargs'].append(dict(color=colors[i], label = g))
                s_arguments['fill_kwargs'].append(dict(color=colors[i], alpha = 0.3))
            
            s_fig = metaseq.plotutils.imshow(array, **s_arguments)
            s_fig = format_default_axes(s_fig, sample = label + ' (' + str(array.shape[0]) + ' features)', feature_name = "Metagenes")
            s_fig = format_subsets_axes(s_fig, by = s_arguments['subset_by'], order = s_arguments['subset_order'], ncols = len(arguments['fill_kwargs']), ticks = ticks, tick_labels = tick_labels, alignment = alignment)

            pp.savefig(s_fig)
            plt.close(s_fig)

            plot_subsets_separate(array= array, arguments = s_arguments, pp = pp, sample = label + ' (' + str(narray.shape[0]) + ' features)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)

            if 'order' in list(fsets.columns.values) and False:# {{{
            # ordering within the cluster is done from bottom to top of heatmap e.g
            # e.g if ordering is the df on the left, on the heatmap from top to bottom will be
            # a 1 classA    classB c
            # b 2 classA    classB d
            # c 3 classB    classA a
            # d 4 classB    classA b
                arguments['sort_by'] = fsets.order
                o_fig = metaseq.plotutils.imshow(array, **o_arguments)
                o_fig = format_default_axes(o_fig, sample = label + ' (' + str(array.shape[0]) + ' features - ordered)', feature_name = "Metagenes", ticks = ticks, tick_labels = tick_labels, alignment = alignment)
                o_fig = format_subseto_axes(o_fig, by = o_arguments['subset_by'], order = o_arguments['subset_order'])
                pp.savefig(o_fig)
                plt.close(o_fig)# }}}

        return {'narray': array, 'features': features, 'fsets': fsets, 'fset_order':fset_order}
        # }}}

# }}}

# For each bam calculate signal and plot it# {{{
# The windows we'll get signal over
def main():
    normalised_arrays = OrderedDict()
    # Create genomic_signal objects that point to data files
    for f in files:
        npz = os.path.splitext(f)[0]
        # Load the windows and arrays
        from metaseq import persistence
        features, arrays = persistence.load_features_and_arrays(prefix = npz)
       
        if 'bw' in arrays:
            normalised_arrays[os.path.basename(npz)] = arrays['bw']
        else:
            # Normalize IP to the control
            normalised = arrays['ip'] - arrays['input']
            normalised_arrays[os.path.basename(npz)] = normalised
    
    narray = np.concatenate((normalised_arrays[normalised_arrays.keys()[0]], normalised_arrays[normalised_arrays.keys()[1]], normalised_arrays[normalised_arrays.keys()[2]]), axis=1)

    from metaseq.results_table import ResultsTable
    features = os.path.splitext(files[1])[0] + '.features'
    features = pybedtools.BedTool(features)

    fsets = args.subsets
    fsets = ResultsTable(fsets, import_kwargs=dict(index_col=0))
    label = os.path.splitext(normalised_arrays.keys()[0])[0]
    test = plot_metagene(narray = narray, features = features, label = label, fsets = fsets)

# }}}

main()
