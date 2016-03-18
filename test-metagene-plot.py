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

files=[ '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.genes-filtered.npz',
        '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.5000-gene_start.npz',
        '/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/mm10/profiles/7E12_2i-M2_0_ddup.gene_end-5000.npz']

# Plotting metagene# {{{
def plot_metagene(narray, subsets = None, features = None, expression = None, label = None, test = None):
    
    def format_axes(figure, name, subsets = False, by= None, order = None, expression = None, test = None):# {{{
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

        if expression is not None:
            print "Formatting expression panel"
            if len(expression) == len(narray):
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
                
                if False:# {{{
                    # log1p = log(1+x)  baseMean = exp(log1p) -1
                    log_expression = np.log1p(expression.baseMean)

                    # Calculating quantiles
                    n_quantiles = 4
                    # True/False array, False are NaN values for baseMean or 0
                    # if I don't remove 0 then when making the qcut I have non-unique bins
                    valid = (np.isfinite(expression.baseMean) & expression.baseMean > 0).values
                    import pandas
                    # Discretize variable into equal-sized buckets based on rank or based on sample quantiles.
                    # http://stackoverflow.com/questions/32011359/convert-categorical-data-in-pandas-dataframe
                    qlabels = pandas.qcut(log_expression[valid], n_quantiles).cat.codes
                    quantiles = (qlabels + 1) / float(n_quantiles)
                    expression.data['expression_quantile'] = 0
                    expression.data['expression_quantile'].loc[valid] = quantiles
                    qlabel = {
                                0: 'silent',
                                0.25: '0-25%',
                                0.5: '25-50%',
                                0.75: '50-75%',
                                1: '75-100%'}

                    # Create linearly spaced colors along a colormap, but with some padding at
                    # the beginning to avoid too-light colors.
                    quantile_colors = matplotlib.cm.YlOrBr(np.linspace(.5, 1, n_quantiles + 1))
                    uquantiles = sorted(list(set(quantiles)))
                    # Plot expression data in the right-most panel
                    figure.strip_axes.fill_betweenx(
                            x1 = np.log1p(sorted(expression.baseMean)),
                            #x1 = sorted(expression.baseMean),
                            y = np.arange(len(expression)),
                            x2=0,
                            color='.5')
                        # Plot the quantile lines as different colors
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
                    # }}}
            else:
                print "FPKMs for %s genes. Genes in array are %s. Not printing expression panel." % (len(expression), len(narray))

        if subsets is not None:
            if by is not None and order is not None:
                metaseq.plotutils.add_labels_to_subsets(
                        figure.array_axes,
                        subset_by = by,
                        subset_order = order)
                figure.line_axes.legend(loc='best', frameon=False)

                # Create another axes along the bottom.
                bottom_ax = plt.subplot(figure.gs[2, 0])
    

        figure.subplots_adjust(left=0.2, right=0.8, hspace=0.05, top=0.95)

        # Positioning the colobar axes took some tweaking:
        figure.cax.set_position([0.75, .15, .05, .10])
        figure.set_facecolor('w')

        return(figure)
# }}}

    arguments = {
            'vmin':5, 'vmax':95, 'percentile':True,
            'figsize':(5, 8),
            # tweak the line plot
            'line_kwargs':dict(color = 'k'),
            'fill_kwargs':dict(color='k', alpha=0.4)
            }

    plt.rcParams['font.family'] = 'Sans'
    plt.rcParams['font.size'] = 10

    window = np.linspace(-5000, 0, 500)
    window = np.append(window, np.linspace(0,100,100))
    window = np.append(window, np.linspace(100,5000,500))
    arguments['x'] = window

    if expression is not None:
        # log1p = log(1+x)  baseMean = exp(log1p) -1
        log_expression = np.log1p(expression.baseMean)
        arguments['sort_by'] = log_expression

        # Params for the strip plot alongside the array
        arguments['subplot_params'] = dict(wspace=0.1, hspace=0.1, bottom=0.05),
        arguments['strip'] = True
        arguments['width_ratios'] = (4, 1)
        arguments['height_ratios'] = (4, 1, 1)
    else:
        arguments['sort_by'] = narray.mean(axis=1)

    pp = PdfPages('/nfs/research2/bertone/user/mxenoph/hendrich/chip/hendrich_2013/' + label + '-profiles.pdf')
    fig = metaseq.plotutils.imshow(narray, **arguments)
    print expression
    fig = format_axes(fig, name = label, expression = expression)

    pp.savefig(fig)
    plt.close(fig)
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
expression = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/inducible/results/report/expression_test_fake_fpkm.tsv'
expression = ResultsTable(expression, import_kwargs=dict(index_col=0))
print "expression len:%s features len:  %s" % (len(expression), len(features))
expression = expression.reindex_to(features, attribute='gene_id')
print "After reindexing: expression len:%s features len:  %s" % (len(expression), len(features))

narray = np.concatenate((normalised_arrays['7E12_2i-M2_0_ddup.5000-gene_start'], normalised_arrays['7E12_2i-M2_0_ddup.genes-filtered'], normalised_arrays['7E12_2i-M2_0_ddup.gene_end-5000']), axis=1)
label = os.path.splitext(normalised_arrays.keys()[0])[0]
#plot_metagene(narray = narray, expression=expression, label = label )
plot_metagene(narray = narray, label = label, test = expression)
