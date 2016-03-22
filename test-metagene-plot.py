#!/usr//bin/env/ python

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
from gffutils.helpers import asinterval# }}}

warnings.simplefilter("error")

files=[ '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.genes-filtered.npz',
        '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.5000-gene_start.npz',
        '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.gene_end-5000.npz']

def get_N_HexCol(N=15):# {{{
    import colorsys
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#' + "".join(map(lambda x: chr(x).encode('hex'),rgb)))
    return hex_out
# }}}

# Plotting metagene# {{{
def plot_metagene(narray, features, label = None, subsets = None):
    def format_axes(figure, name, by = None, order = None, expression = None, quantiles = None):# {{{
        figure.line_axes.axis('tight')
        for ax in [figure.line_axes, figure.array_axes]:
            ax.axvline(0, color='k', linestyle='--')
            
        # Label axes
        figure.array_axes.set_title(name)
        figure.array_axes.set_ylabel("Genes")

        figure.line_axes.set_ylabel("Average enrichment (IP - input)\n RPMMR")
        figure.line_axes.set_xlabel("Metagenes ")

        # Line at the TTS -- grey
        figure.line_axes.axvline(0, color='k', linestyle=':')
        figure.line_axes.axvline(100, color='c', linestyle=':')

        figure.cax.set_ylabel("RPMMR")
        
        # Remove ticks from heatmap x axis
        for tick in figure.array_axes.get_xticklabels():
            tick.set_visible(False)

        # Add average plot per subsets
        if by is not None and order is not None:# {{{
            print "In format_axes subsets"
            metaseq.plotutils.add_labels_to_subsets(
                    figure.array_axes,
                    subset_by = by,
                    subset_order = order)
            figure.line_axes.legend(loc ='best', frameon = False)
        # }}}

        if expression is not None:# {{{
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

            figure.strip_axes.set_xlabel('log(FPKM)')
            figure.strip_axes.xaxis.set_major_locator(
                matplotlib.ticker.MaxNLocator(nbins = 2, prune = None))

            # Plot expression data in the right-most panel
            figure.strip_axes.fill_betweenx(
                    x1 = sorted(expression),
                    y = np.arange(len(expression)),
                    x2 = 0,
                    color = '.5')

            # If information about silent, Q1-Q4 genes is given add a subplot to the figure
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
                #    bottom_ax.set_xticks(ticks)
                #    bottom_ax.set_xticklabels(ticklabels)
                #    for txt, ha in zip(bottom_ax.get_xticklabels(), alignment):
                #        txt.set_horizontalalignment(ha)
                #    bottom_ax.axis('tight')
                #    bottom_ax.axvline(100, **border_style)
                #    bottom_ax.axvline(600, **border_style)
                #    bottom_ax.legend(loc='best', fontsize=8, frameon=False, labelspacing=0.4)
                #    bottom_ax.set_ylabel('Enrichment\n(RPMMR)')# }}}

        figure.subplots_adjust(left = 0.2, right = 0.8, hspace = 0.05, top = 0.95)

        # Positioning the colobar axes took some tweaking:
        figure.cax.set_position([0.75, .15, .05, .10])
        figure.set_facecolor('w')

        return(figure)
# }}}

    arguments = {
            'vmin':5, 'vmax':95, 'percentile':True,
            'figsize':(5, 8)
            }

    plt.rcParams['font.family'] = 'Sans'
    plt.rcParams['font.size'] = 10

    window = np.linspace(-5000, 0, 500)
    window = np.append(window, np.linspace(0,100,100))
    window = np.append(window, np.linspace(100,5000,500))
    arguments['x'] = window

    pp = PdfPages('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/' + label + '-profiles.pdf')

    if subsets is not None:
        assert type(subsets) is not 'metaseq.results_table.ResultsTable', \
                'Subsets passed is not a table.'

        assert len([i for i in ['groups', 'expression', 'order'] if i in list(subsets.columns.values)]) != 0, \
                'Table does not contain any column that I know how to handle'

        subsets_original_order = subsets.index
        if len(features) != len(subsets): # {{{
            tmp = [g.name for g in features if g.name not in subsets.index]
            print "%s features in array are not present in subseting file. These will be presented as the UNK group when grouping or ignored when ordering heatmap." % len(tmp)
            # Joining the 2 tuples. tmp contains the genes to be filtered out
            tmp = [x for x in subsets.index if x not in [g.name for g in features]]
            print "%s features present in subseting file but not array. These will not be plotted." % len(tmp)
            
            features_ids = [g.name for g in features]
            # Filtering based on original order as reindexing in line 173 will result in all features in array been represented in subsets
            array_indices_to_keep = [i for i, gene in enumerate(features) if gene.name in subsets_original_order]
            narray_subset = narray[array_indices_to_keep]
            features_subset = features.filter(lambda gene: gene.name in subsets.index).saveas()
            subsets_set = subsets.reindex_to(features_subset, attribute = 'gene_id')
            print "subsets %s subsets_set %s" % (len(subsets), len(subsets_set))

            # zip(a, b) combines two equal-length lists and meregs tehm together in pairs. e.g a = ['a', 'b'] and b = [1, 2] zip(a,b) = [('a',1), ('b',2)]
            # will return True only for those indeces where the gene_id in both features and subset file are the same. If subsets is ordered as features 
            # will return True for len(features) = len(subsets)
            if len(features) != len([True for i,j in zip(subsets.index, [gene.name for gene in features]) if i==j]):
                print "Features provided in subset file are not ordered as array. Re-indexing"
                subsets = subsets.reindex_to(features, attribute = 'gene_id')
        # }}}

        if 'groups' in list(subsets.columns.values) and False: # {{{
            # Subsetting based on groups provided
            print "plot_metagene(): subsetting based on groups provided"
            # For getting the groups keep only those indeces that were in the original subset order, otherwise
            # gene_ids in the array but not the subset file will produce nan. That makes it trickier to remove and fix
            # downstream, expecially if I ever want a group to be named nan in the subsets file
            groups = subsets.iloc[subsets.index.isin(subsets_original_order)].groups.unique()
            cls = np.zeros(len(narray)).astype('str')
            
            # The numpy arrays for each group label will have the same length as 
            # len(features) == len(subsets) and same order because before getting here subsets was reindexed
            subset = []
            for g in groups:
                subset.append((g, (subsets.groups == g).values))
            
            if not all(subsets.index.isin(subsets_original_order)):
                print "Features in the array are missing from the subsets file. group those together under the UNK group"
                subset.append(('UNK group', ~(subsets.index.isin(subsets_original_order))))
                groups = np.append(groups, 'UNK group')

            subset = tuple(subset)
            for lab, ind in subset:
                cls[ind] = lab
            # if group names are numbers then if 0 is a group asset will fail even if all features are assigned to a group
            # Commenting this out until I can think of a better way to do it
#            assert sum(cls == '0.0') == 0

            arguments_subset  = arguments
            # Saving groups and subsets to arguments for plotutils
            arguments_subset['subset_by'] = cls
            arguments_subset['subset_order'] = sorted(groups)

            arguments_subset['line_kwargs'] = []
            arguments_subset['fill_kwargs'] = []
            colors = get_N_HexCol(len(groups))
            
            for i, g in enumerate(groups):
                arguments_subset['line_kwargs'].append(dict(color=colors[i], label = g))
                arguments_subset['fill_kwargs'].append(dict(color=colors[i], alpha = 0.3))
            
            arguments_subset['sort_by'] = narray.mean(axis=1)
            fig_subset = metaseq.plotutils.imshow(narray, **arguments_subset)
            fig_subset = format_axes(fig_subset, name = label, by = arguments_subset['subset_by'], order = arguments_subset['subset_order'])
            pp.savefig(fig_subset)
            plt.close(fig_subset)
        # }}}
        
        if 'expression' in list(subsets.columns.values) and False:# {{{
            print "Printing expression strip"
            arguments_expression = arguments
            arguments_expression['sort_by'] = narray.mean(axis=1)

            # Params for the strip plot alongside the array
            arguments_expression['subplot_params'] = dict(wspace=0.1, hspace=0.1, bottom=0.05)
            arguments_expression['strip'] = True
            arguments_expression['width_ratios'] = (4, 1)
            arguments_expression['height_ratios'] = (4, 1, 1)
            
            fig_expression = metaseq.plotutils.imshow(narray, **arguments_expression)
            #fig_expression = format_axes(fig_expression, name = label, expression = subsets.expression, quantiles = subsets.quartiles_valid)
            fig_expression = format_axes(fig_expression, name = label, expression = subsets.expression)
            pp.savefig(fig_expression)
            plt.close(fig_expression)
        # }}}
        if 'order' in list(subsets.columns.values):
            print "Order column provided...Subsetting and/or ordering"
#            features_ids = [g.name for g in features]
#            # Filtering based on original order as reindexing in line 173 will result in all features in array been represented in subsets
#            array_indices_to_keep = [i for i, gene in enumerate(features) if gene.name in subsets_original_order]
#            narray_subset = narray[array_indices_to_keep]
#            features_subset = features.filter(lambda gene: gene.name in subsets.index).saveas()
#            subsets_set = subsets.reindex_to(features_subset, attribute = 'gene_id')
#            print "subsets %s subsets_set %s" % (len(subsets), len(subsets_set))
#
            arguments['line_kwargs'] = dict(color = 'k')
            arguments['fill_kwargs'] = dict(color='k', alpha=0.4)
            arguments['sort_by'] = subsets_set.order
#            arguments['sort_by'] = narray_subset.mean(axis=1)
            fig = metaseq.plotutils.imshow(narray_subset, **arguments)
            fig = format_axes(fig, name = label)
            pp.savefig(fig)
            plt.close(fig)
    else:
        print "No csv file provided for subseting or ordering"
        arguments['line_kwargs'] = dict(color = 'k')
        arguments['fill_kwargs'] = dict(color='k', alpha=0.4)
        fig = metaseq.plotutils.imshow(narray, features = features_subset, **arguments)
        fig = format_axes(fig, name = label)
        pp.savefig(fig)
        plt.close(fig)
    # }}}

    pp.close()
# }}}

from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict
normalised_arrays = OrderedDict()

for f in files:
    npz = os.path.splitext(f)[0]
    # Load the windows and arrays
    from metaseq import persistence
    # Only keeping the last BedTool object for features but since the order is the same in all 3 that is ok
    features, arrays = persistence.load_features_and_arrays(prefix = npz)
   
    # Normalize IP to the control
    normalised = arrays['ip'] - arrays['input']
    normalised_arrays[os.path.basename(npz)] = normalised
    
from metaseq.results_table import ResultsTable
# if I pass the gtf then attribute gene_id can not be found. weird but that's how it is. Should Potentially change all suffices to gff
features = pybedtools.BedTool('/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.metaseq.genes.filtered.gtf')
expression = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/results/report/expression_test_fake_fpkm2.tsv'
#expression = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/results/report/ordered_subset.tsv'
expression = ResultsTable(expression, import_kwargs=dict(index_col=0))
print "Subset file len(features): %s len(annotation):  %s" % (len(expression), len(features))
#expression = expression.reindex_to(features, attribute='gene_id')
#print "After reindexing, len(features): %s len(annotation):  %s" % (len(expression), len(features))

narray = np.concatenate((normalised_arrays['7E12_2i-M2_0_ddup.5000-gene_start'], normalised_arrays['7E12_2i-M2_0_ddup.genes-filtered'], normalised_arrays['7E12_2i-M2_0_ddup.gene_end-5000']), axis=1)
label = os.path.splitext(normalised_arrays.keys()[0])[0]
#plot_metagene(narray = narray, expression=expression, label = label )
plot_metagene(narray = narray, features = features, label = label, subsets = expression)
