#!/usr//bin/env/ python

# Import modules# {{{
# Needed for division to return float
from __future__ import division
import os, sys, argparse, re
from time import gmtime, strftime
import fnmatch
import gffutils
import numpy as np
import pybedtools
import metaseq
from pylab import *
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval# }}}

files=[ '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.5000-gene_start.npz',
        '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.genes-filtered.npz',
        '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.gene_end-5000.npz']
expression = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/results/report/expression_test_fake_fpkm.tsv'

# Plotting metagene# {{{
def plot_metagene(narray, subsets = None, features = None, feature_type = "Genes", expression = None):
    
    def format_axes(figure, name, subsets = False, by= None, order = None, expression = None):# {{{
        figure.line_axes.axis('tight')
        for ax in [figure.line_axes, figure.array_axes]:
            ax.axvline(0, color='k', linestyle='--')
            
        # Label axes
        figure.array_axes.set_title(name)
        figure.array_axes.set_ylabel("Genes")
        figure.line_axes.set_ylabel("Average enrichment (IP - input)\n RPMMR")
        figure.line_axes.set_xlabel("Metagenes ")
        figure.cax.set_ylabel("Enrichment (IP - input)\n RPMMR")

        if subsets:
            if by is not None and order is not None:
                metaseq.plotutils.add_labels_to_subsets(
                        figure.array_axes,
                        subset_by = by,
                        subset_order = order)
                figure.line_axes.legend(loc='best', frameon=False)
                
        # Annotate the array axes
#        ticks, ticklabels, alignment = zip(*[
#            (0, '5kb\nupstream', 'right'),
#            (100, 'TSS', 'center'),
#            (350, 'gene body', 'center'),
#            (600, 'polyA\nsite', 'center'),
#            (700, '5kb \ndownstream', 'left'),
#        ])
#        figure.line_axes.set_xticks(ticks)
#        for txt in figure.array_axes.get_xticklabels():
#            txt.set_visible(False)

#
    # Axes limit tweaks
#    figure.strip_axes.axis('tight')
#    figure.strip_axes.axis(xmin=0)
#
#    # Remove all edges but the bottom x-axis.
#    figure.strip_axes.yaxis.set_visible(False)
#    figure.strip_axes.spines['top'].set_visible(False)
#    figure.strip_axes.spines['right'].set_visible(False)
#    figure.strip_axes.spines['left'].set_visible(False)
#    figure.strip_axes.xaxis.set_ticks_position('bottom')
#
#    figure.strip_axes.set_xlabel('log(RPKM)')
#    figure.strip_axes.xaxis.set_major_locator(
#        matplotlib.ticker.MaxNLocator(nbins=4, prune=None))
#
#    # First line axes --------------------------------------------------------
#    figure.line_axes.set_ylabel('Enrichment\n(RPMMR)')
#    figure.line_axes.axvline(100, **border_style)
#    figure.line_axes.axvline(600, **border_style)
#    figure.line_axes.set_xticklabels([])
#    figure.line_axes.set_xticks(ticks)
#
#    figure.subplots_adjust(left=0.2, right=0.8, hspace=0.05, top=0.95)
#
#    # Annotate and resize the colorbar axes
#    figure.cax.set_ylabel('Enrichment\n(RPMMR)')
#
#    # Positioning the colobar axes took some tweaking:
#    figure.cax.set_position([0.75, .15, .05, .12])
#
#    # Create another axes along the bottom.
#    bottom_ax = plt.subplot(figure.gs[2, 0])
#
#    # Plot the quantile lines as different colors
#    for q, color in zip(uquantiles, quantile_colors)[::-1]:
#        ind = (e.data['expression_quantile'] == q).values
#        plotutils.ci_plot(
#            np.arange(700),
#            diffed[ind, :],
#            ax=bottom_ax,
#            line_kwargs=dict(
#                color=color, label=qlabel[q], alpha=0.5),
#            fill_kwargs=dict(
#                color=color, alpha=0.3),
#        )
#
#    # Axes configure on the bottom axes.
#    bottom_ax.set_xticks(ticks)
#    bottom_ax.set_xticklabels(ticklabels)
#    for txt, ha in zip(bottom_ax.get_xticklabels(), alignment):
#        txt.set_horizontalalignment(ha)
#    bottom_ax.axis('tight')
#    bottom_ax.axvline(100, **border_style)
#    bottom_ax.axvline(600, **border_style)
#    bottom_ax.legend(loc='best', fontsize=8, frameon=False, labelspacing=0.4)
#    bottom_ax.set_ylabel('Enrichment\n(RPMMR)')
#
#    # Final figureure setup, add the subplot labels, and save.
#    figure.set_facecolor('w')
#    filename = os.path.join(figuredir, ip_label + '.png')
#    subplot_label_kwargs = dict(size=14, weight='bold',
#                                verticalalignment='top')
#    for ax, letter, pos in [
#        (figure.array_axes, 'A', (-.3, 1.)),
#        (figure.line_axes, 'B', (-.3, 1.)),
#        (bottom_ax, 'C', (-.3, 1.)),
#    ]:
#        ax.text(pos[0], pos[1], letter, transform=ax.transAxes,size=14,
#                weight='bold', verticalalignment='top')
#

        return(figure)# }}}

    arguments = {
            'vmin':5, 'vmax':95, 'percentile':True,
            'figsize':(5, 8),
            }

    if features is not None and subsets is not None:# {{{
        print "Subsetting based on file %s" % subsets
        from metaseq.results_table import ResultsTable
        subsets = ResultsTable(subsets, import_kwargs=dict(index_col=0))
        # file needs to have 2 columns, gene_id and group
        groups = subsets.groups.unique()
        features_from_groups = subsets.index
        subsets = subsets.reindex_to(features, attribute='gene_id')
        cls = np.zeros(len(array)).astype('str')

        subset = []
        for g in groups:
            subset.append((g, (subsets.groups == g).values))
        
        subset.append(('UNK group', ~(subsets.index.isin(features_from_groups))))
        subset = tuple(subset)

        for label, ind in subset:
            cls[ind] = label
        assert sum(cls == '0.0') == 0

        # Features found in the gtf used to construct the array 
        # but not described in subsets are set to UNK group
        groups = np.append(groups, 'UNK group')
        # Saving groups and subsets to arguments for plotutils
        arguments['subset_by'] = cls
        arguments['subset_order'] = sorted(groups)

        arguments['line_kwargs'] = []
        arguments['fill_kwargs'] = []
        colors = get_N_HexCol(len(groups))

        for i, g in enumerate(groups):
            arguments['line_kwargs'].append(dict(color=colors[i], label = g))
            arguments['fill_kwargs'].append(dict(color=colors[i], alpha = 0.3))

        #arguments['sort_by'] = subsets.index
        # gs = gene_start.filter(lambda b: b.name in sets.index)
        # print "N of TSS in subset:", len(gs)
        # }}}
    
    if features is not None and expression is not None: # {{{
        print "Subsetting based on expression (%s)" % expression
        deseq_results = metaseq.results_table.DESeq2Results(expression)
        deseq_results = deseq_results.reindex_to(features, attribute = 'gene_id')
        cls = np.zeros(len(array)).astype('str')

#        print 'Type: %s, length of features: %s' % (feature_type, len(array))
        # This is a tuple (just like list though can't be changed and defined with parenthesis)
        # each ind is labeled. All 3 ind have the same length, equal to the number of genes in deseq_results
        # At this point deseq_results and features used for numpy array should have same length and order
        subset = (
                ('unchanged', deseq_results.unchanged(0.05).values),
                ('down', deseq_results.downregulated(0.05).values),
                ('up', deseq_results.upregulated(0.05).values))

        for label, ind in subset:
            cls[ind] = label

        # Make sure all genes are classified in the 3 categories and none is left 0 from
        # the initial declaration
        assert sum(cls == '0.0') == 0

        # Saving groups and order to arguments for plotutils
        arguments['subset_by'] = cls
        arguments['subset_order'] = ['unchanged', 'down', 'up']

        arguments['line_kwargs'] = [dict(color='#f57900', label = 'up'), 
                            dict(color='#8f5902', label = 'down'),
                            dict(color='#000000', label = 'unchanged')]
        arguments['fill_kwargs'] = [dict(color='#f57900', alpha = 0.3),
                            dict(color='#8f5902', alpha = 0.3),
                            dict(color='#000000', alpha = 0.3)]
    elif subsets is None:
        print 'Expression file not provided'
        arguments['line_kwargs'] = dict(color='k')
        arguments['fill_kwargs'] = dict(color='k', alpha=0.3)
        # }}}

#    if feature_type == 'genes':
#        arguments['imshow_kwargs'] = dict(interpolation = 'none')

    # narray is the normalised_arrays dictionary
    blocks = ['5000-gene_start', 'genes', 'gene_end-5000']
    
    names = set()
    for key in narray.keys():
        print key
        basename = os.path.splitext(os.path.basename(key))[0]
        names.add(basename)

    for n in names:
        print "NAME", n
        for key in narray.keys():
            print key
            print narray.get(key).shape
        print [n + '.' + s for s in blocks]
        #metagene = np.column_stack(map(narray.get, [n + '.' + s for s in blocks]))
        test  = map(narray.get, [n + '.' + s for s in blocks])

        # Plot expression data in the right-most panel
    if expression is not None:
        assert len(expression) == len(metagene), \
                "FPKMs for %s genes. Genes in array are %s" % (len(expression), len(metagene))

        figure.strip_axes.fill_betweenx(
                x1 = np.log1p(sorted(rpkm)),
                y = np.arange(len(rpkm)),
                x2=0,
                color='.5')

        window = np.linspace(-2000, 500, 100)
        window = np.append(window, np.linspace(0,100,100))
        arguments['x'] = window
        arguments['sort_by'] = metagene.mean(axis=1)
        arguments['strip'] = True
        arguments['subplot_params'] =dict(wspace=0.1, hspace=0.1, bottom=0.05)
        arguments['width_ratios'] = (4, 1)
        arguments['height_ratios'] = (4, 1, 1)

        pp = PdfPages('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/' + n + '-profiles.pdf')
        fig = metaseq.plotutils.imshow(metagene, **arguments)
        fig = format_axes(fig, name = n, expression =expression)

        pp.savefig(fig)
        plt.close(fig)
        pp.close()
        return(metagene)
        # }}}

from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict
normalised_arrays = OrderedDict()
features = dict()

# Files should be ordered upstream, genes, downstream for np_column_stack to work as intended
# Using OrderedDict to preserve the order of the passed arrays
for f in files:
    npz = os.path.splitext(f)[0]
    # Load the windows and arrays
    from metaseq import persistence
    features, arrays = persistence.load_features_and_arrays(prefix = npz)
   
    # Normalize IP to the control
    normalised = arrays['ip'] - arrays['input']
    normalised_arrays[os.path.basename(npz)] = normalised
    
from metaseq.results_table import ResultsTable
expression = ResultsTable(expression, import_kwargs=dict(index_col=0))
narray = np.column_stack(normalised_arrays)
#expression = expression.reindex_to(features, attribute='gene_id')
plot_metagene(narray = narray, expression=expression)
