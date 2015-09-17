#!/usr/bin/env python

# ---------------------------------------------------------
# Get TSS +/- 1kb for all annotated transcripts.
#

# Import modules# {{{
# Needed for division to return float
from __future__ import division
import os, sys, argparse, re
from time import gmtime, strftime
import gffutils
import numpy as np
import pybedtools
import metaseq
from pylab import *
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval# }}}

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()
# Required arguments
parser.add_argument('-p', '--path', metavar="path", type=str, required=True, help = 'Path to search for npz files')
# Default arguments
# Optional arguments
#parser.add_argument('-b', '--bound', metavar="file", type=str, required=False, help= 'TSV file with gene_id (ensembl) and bound (1/0) fields')
#parser.add_argument('-s', '--order', type=str, required=False, help= 'Keep order of features provided')

args = parser.parse_args()

if not args.path.endswith(os.sep):
    args.path = args.path + os.sep

plot_path = args.path + 'plots/'

# Create plot directory if it doesn't exist
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

files = [ os.path.join(args.path, f) for f in os.listdir(args.path) if re.match(r'.*\.npz',f) ]

#}}}

# Functions# {{{
class ref:
    def __init__(self, obj): self.obj = obj
    def get(self):    return self.obj
    def set(self, obj):      self.obj = obj

def get_N_HexCol(N=15):# {{{
    import colorsys
    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    hex_out = []
    for rgb in HSV_tuples:
        rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        hex_out.append('#' + "".join(map(lambda x: chr(x).encode('hex'),rgb)))
    return hex_out# }}}

def most_common(lst):
    return max(set(lst), key=lst.count)

# Plotting average enrichment -- geom_line# {{{
def plot_average(array, window, pp, name = None, feature_type = 'tss'):
    from matplotlib import pyplot as plt
    plt.rcParams['font.size'] = 11
    plt.rcParams['legend.scatterpoints'] = 1
    plt.rcParams['legend.fontsize'] = 10
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # for plotting more than one IP on the same plot, array should be a dictionaty
    # with keys the basename and values the normalised array
    if isinstance(array, dict):
        n_samples = len(array.keys())
        colors = get_N_HexCol(n_samples)

        for key in array.keys():
            ax.plot(window, array.get(key).mean(axis=0), color=colors[array.keys().index(key)], label = key)
    elif name is not None:
        ax.plot(window, array.mean(axis=0), color = 'r', label = name)
        # Add a vertical line at TSS/gene-start
        ax.axvline(0, color='k', linestyle='--')
    else:
        print "You did not pass a name for the sample."
    
    # Add labels and legend
    ax.set_xlabel('Distance from ' + feature_type +' (bp)')
    ax.set_ylabel('Normalised average read coverage (per M mapped reads)')
    ax.legend(loc = 2, frameon = False, fontsize = 14, labelspacing = .3, handletextpad = 0.2)
    pp.savefig(fig)
    # }}}

# bound should be a all genes with 1 indicating bound and 0 not bound
# TODO: replicates should be plotted together split by bound and not bound
# add check if expression matrix provided and if so call signal_expr
# kwargs = optional arguments
def signal_bound(norm_sub, subset, window, pp, **kwargs):# {{{
   # If expr object/dataframe not provided then set to none for the next checks
    expr = kwargs.get('expr', None)

    import numpy as np
    from matplotlib import pyplot as plt

    window = window / 2
    #x = np.linspace(-window, window, 100)
    x = np.linspace(-2000, 1000, 100)
    
    n_samples = len(norm_sub.keys())
    colors = get_N_HexCol(n_samples)

#    for key in norm_sub.keys():
#        normalized_subtracted = norm_sub.get(key)
#        
#        fig_simple = metaseq.plotutils.imshow(
#                normalized_subtracted,
#                x=gene_start,
#                # Increase the contrast by truncating the colormap
#                # to 5th and 95th percentiles
#                vmin=5,
#                vmax=95,
#                percentile=True,
#                sort_by=normalized_subtracted.mean(axis=1),
#                # tweak the line plot
#                line_kwargs=dict(color='k'),
#                fill_kwargs=dict(color='k', alpha=0.4))
#        
#        # auto-zoom axes
#        fig_simple.line_axes.axis('tight')
#        
#        # draw a vertical line at zero on both axes
#        for ax in [fig_simple.line_axes, fig_simple.array_axes]:
#            ax.axvline(0, color='k', linestyle='--')
#            # Label axe
#            fig_simple.array_axes.set_ylabel("Genes")
#            fig_simple.line_axes.set_ylabel("Average enrichment\n(IP - input)\nreads per million mapped reads")
#            fig_simple.line_axes.set_xlabel("Distance from TSS (bp)")
#            fig_simple.cax.set_ylabel("Enrichment\n(IP - input)\nreads per million mapped reads)");
#            pp.savefig(fig_simple)
#
    fig = plt.figure()
    # For all bound TSS
    ax = fig.add_subplot(111)
    heatmaps = dict()
    norm_sub_bound = dict()
    
    gs = gene_start.filter(lambda: b.name in set.index)
    # Converting numbers to strings and making a matrix to hold the information for bound/not# {{{
    subset_by = [str(num) for num in subset.clusters]
    subset_by = array(subset_by, dtype="|S32")

    subset_order = ['3', '2', '1', '4', '5']
    # }}}

    for key in norm_sub.keys():
        normalized_subtracted = norm_sub.get(key)

        heatmaps[key]= metaseq.plotutils.imshow(normalized_subtracted,
                x=x,
                features=gs,
                vmin=5, vmax=95, percentile=True,
                # axis=1 means that mean is calculated on the rows
#                sort_by=normalized_subtracted.mean(axis=1),
                subset_by= subset_by,
                subset_order=subset_order,
                line_kwargs=[dict(color='k', label='1'),
                    dict(color='r', label='2'),
                    dict(color='r', label='3'),
                    dict(color='r', label='4'),
                    dict(color='r', label='5')],
                fill_kwargs=[dict(color='k', alpha=0.3),
                    dict(color='r', alpha=0.3),
                    dict(color='r', alpha=0.3),
                    dict(color='r', alpha=0.3),
                    dict(color='r', alpha=0.3),
                    ]
                )
        # adding labels to heatmap
        metaseq.plotutils.add_labels_to_subsets(heatmaps[key].array_axes, subset_by=subset_by, subset_order=subset_order,)
        # adding labels to average signal plot
        heatmaps[key].line_axes.legend(loc='best', frameon=False)
#        heatmaps[key].suptitle(key + "\n(Bound genes: " + str(len(subset[subset_by=='Bound'])) + ' )')
        pp.savefig(heatmaps[key])

        # bound.bound ==1 returns a TRUE/FALSE vector and it's in the same order as the tsses
        # that used to make the normalized_subtracted array, hence order is the same
        normalized_subtracted = normalized_subtracted[(subset.boolean==1).values, :]
        color = colors[norm_sub.keys().index(key)]

        metaseq.plotutils.ci_plot(x,
                normalized_subtracted,
                ax = ax,
                line_kwargs=dict(color= color, label=key),
                fill_kwargs=dict(color= color, alpha=0.3))

        norm_sub_bound[key]= normalized_subtracted

    ax.set_xlabel('Distance from ' + 'TSS' +' (bp)')
    ax.set_ylabel('Normalised average read coverage (per million mapped reads)')
    ax.set_title('Subset: Bound')
    ax.legend(loc=2, frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
    fig.tight_layout()
    pp.savefig(fig)

    try:
        expr
    # if expr is None that means no argument expr was provided to this function
    except None:
        print "Plotting signal only for bound genes."
    else:
        expr_subset = expr.loc[bound[subset.boolean==1].index]
        axes = signal_expr(norm_sub_bound, expr_subset, window*2, ' (Bound + DE)')
        pp.savefig(axes);

#    return (fig, axes)


# }}}

def plot_tss(array, window, name, pp, order = None, features = None, feature_type = "Genes", expression = None):# {{{

    def format_axes(figure, feature_type = 'tss', subsets = False, by= None, order = None):# {{{
        figure.line_axes.axis('tight')
        for ax in [figure.line_axes, figure.array_axes]:
            ax.axvline(0, color='k', linestyle='--')
            
        # Label axes
        figure.array_axes.set_ylabel(feature_type)
        figure.line_axes.set_ylabel("Average enrichment (IP - input)\n RPMMR")
        figure.line_axes.set_xlabel("Distance from " + feature_type + " (bp)")
        figure.cax.set_ylabel("Enrichment (IP - input)\n RPMMR")

        if subsets:
            if by is not None and order is not None:
                metaseq.plotutils.add_labels_to_subsets(
                        figure.array_axes,
                        subset_by = by,
                        subset_order = order)
                figure.line_axes.legend(loc='best', frameon=False)
                
        if feature_type == 'genes':
            fig.line_axes.set_xlabel("Percentage of gene")

        return(figure)# }}}

    arguments = {
            'x':window,
            'vmin':5, 'vmax':95, 'percentile':True,
            'figsize':(10, 15),
            }

    if order is None:
        arguments['sort_by'] = array.mean(axis=1)
        print "No order provided"
    else:
        from metaseq.results_table import ResultsTable
        order = ResultsTable(args.order, import_kwargs=dict(index_col=0))
        order = order.reindex_to(features, attribute='gene_id')
        arguments['sort_by'] = order.index
       # gs = gene_start.filter(lambda b: b.name in sets.index)
       # print "N of TSS in subset:", len(gs)
    
    if features is not None and expression is not None: # {{{
        print "Subsetting based on expression"
        deseq_results = metaseq.results_table.DESeq2Results(expression)
        deseq_results = deseq_results.reindex_to(features, attribute = 'gene_id')
        cls = np.zeros(len(array)).astype('str')

        print 'Type: %s, length of features: %s' % (feature_type, len(array))
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
    else:
        print 'Features not provided'
        arguments['line_kwargs'] = dict(color='k')
        arguments['fill_kwargs'] = dict(color='k', alpha=0.3)
        # }}}

    if feature_type == 'genes':
        arguments['vmin'] = 1
        arguments['vmax'] = 99
        arguments['imshow_kwargs'] = dict(interpolation = 'none')

    fig = metaseq.plotutils.imshow(array, **arguments)


    if features is not None and expression is not None:
        fig = format_axes(fig, feature_type = feature_type, subsets = True, by = arguments['subset_by'], order = arguments['subset_order'])
    else:
        fig = format_axes(fig, feature_type = feature_type)


    pp.savefig(fig)# }}}

# }}}

# # For each bam calculate signal and plot it# {{{
# # The windows we'll get signal over
def main():
    normalised_arrays = dict()
    # Create genomic_signal objects that point to data files
    for f in files:
        npz = os.path.splitext(f)[0]
        # Load the windows and arrays
        from metaseq import persistence
        features, arrays = persistence.load_features_and_arrays(prefix = npz)
       
        # Normalize IP to the control
        normalised = arrays['ip'] - arrays['input']
        normalised_arrays[os.path.basename(npz)] = normalised
        
    # Set ensures that the values in list are unique
    # http://stackoverflow.com/questions/12897374/get-unique-values-from-a-list-in-python
    names = set()
    for key in normalised_arrays.keys():
        basename = os.path.splitext(os.path.basename(key))[0]
        names.add(basename)
    print names
    tmp = {'genes':dict(), 'tss':dict(), 'gene_start':dict()}
    windows_tmp = {'genes':dict(), 'tss':dict(), 'gene_start':dict()}
    
    for x in names:
        from matplotlib.backends.backend_pdf import PdfPages
        # Create pdfpage object
        pp = PdfPages(plot_path + x + '-profiles.pdf')

        npz = [ k for k in normalised_arrays.keys() if re.match(x,k) ]
        print npz
        
        for n in npz:
            print 'Plotting for %s' % n
            window = os.path.splitext(n)[1].replace(".", "")
            pattern = re.compile('(\d+)-(\w+)-(\d+)')
            # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
            matching = pattern.match(window)
            
            if matching:
                upstream, feature_type, downstream = matching.groups()
                # bins are of size 100bp in count-tags.py
                x = np.linspace(-int(upstream), int(downstream), 100)
                
                # This should only return one file
                features = [ os.path.join(args.path, f) for f in os.listdir(args.path) if re.match(re.escape(n) + r'.features',f) ][0]
                features = pybedtools.BedTool(features)

#                if feature_type == 'gene_start':
#                    plot_tss(array = normalised_arrays.get(n), name = n, window = x,
#                            pp = pp,
#                            features = features,
#                            feature_type = feature_type,
#                            expression = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/deseq/2i_wt-ko_de.tsv')
#
                if feature_type == 'tss':
                    tmp['tss'][n] = normalised_arrays.get(n)
                elif feature_type == 'gene_start':
                    tmp['gene_start'][n] = normalised_arrays.get(n)

            elif window == 'genes':
                feature_type = 'genes'
                # gene array goes from o) to 100% in bins of 100bp
                x = np.linspace(0, 100, 100)

                features = [ os.path.join(args.path, f) for f in os.listdir(args.path) if re.match(re.escape(n) + r'.features',f) ][0]
                features = pybedtools.BedTool(features)

#                plot_tss(array = normalised_arrays.get(n), name = n, window = x,
#                        pp = pp,
#                        features = features,
#                        feature_type = feature_type,
#                        expression = '/nfs/research2/bertone/user/mxenoph/hendrich/rna/mm10/deseq/2i_wt-ko_de.tsv')

                tmp['genes'][n] = normalised_arrays.get(n)

            else:
                print 'Unknown feature type for %s. Exiting now.' % n
                pp.close()
                sys.exit()
                
        pp.close()

    # TODO: for replicates call plot_average
    if True:
        pp = PdfPages(plot_path + 'all' + '-profiles.pdf')

        for key in tmp.keys():
            # key %in% tss, gene_start, genes
            window = [ w.split('.')[1] for w in tmp.get(key).keys() ]
            print window
            if len(set(window)) == 1:
                print 'all the same'
                pattern = re.compile('(\d+)-(\w+)-(\d+)')
                # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
                matching = pattern.match(window[0])
                
                if matching:
                    upstream, feature_type, downstream = matching.groups()
                    # bins are of size 100bp in count-tags.py
                    x = np.linspace(-int(upstream), int(downstream), 100)
                else:
                    # gene array goes from o) to 100% in bins of 100bp
                    x = np.linspace(0, 100, 100)
                    
                plot_average(tmp[key], x, pp = pp, feature_type = key)
        pp.close()
    sys.exit()
        

# }}}

main()
