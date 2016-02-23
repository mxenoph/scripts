#!/usr/bin/env python

# ---------------------------------------------------------
# Plotting function for gene starts, transcript starts and genes
#

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

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()
# Required arguments
parser.add_argument('-p', '--path', metavar = "path", type = str, required = True, help = 'Path to search for npz files')
# Optional arguments
parser.add_argument('-m', '--matching', type=str, required = False, help = 'Pattern by which to retrieve files from the path. e.g. "*3KO*"',
        # default pattern set to all so I don't have to write an if else condition for defining files
        default = "*")
parser.add_argument('-s', '--subsets', type=str, required = False, help = 'Tab delimited file with gene_id, group columns')
parser.add_argument('-e', '--expression', type=str, required = False, help = 'DESeq output file')

args = parser.parse_args()

if not args.path.endswith(os.sep):
    args.path = args.path + os.sep

plot_path = args.path + 'plots/'

# Create plot directory if it doesn't exist
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

files = [ os.path.join(args.path, f) for f in os.listdir(args.path) if f.endswith(".npz") and fnmatch.fnmatch(f, args.matching)]
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
    def color_variant(hex_color, brightness_offset=1):
        # http://chase-seibert.github.io/blog/2011/07/29/python-calculate-lighterdarker-rgb-colors.html
        """ takes a color like #87c95f and produces a lighter or darker variant """
        if len(hex_color) != 7:
            raise Exception("Passed %s into color_variant(), needs to be in #87c95f format." % hex_color)
        rgb_hex = [hex_color[x:x+2] for x in [1, 3, 5]]
        new_rgb_int = [int(hex_value, 16) + brightness_offset for hex_value in rgb_hex]
        new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int] # make sure new values are between 0 and 255
        # hex() produces "0x88", we want just "88"
        return "#" + "".join([hex(i)[2:] for i in new_rgb_int])

    from matplotlib import pyplot as plt
    plt.rcParams['font.size'] = 11
    plt.rcParams['legend.scatterpoints'] = 1
    plt.rcParams['legend.fontsize'] = 10
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # for plotting more than one IP on the same plot, array should be a dictionaty
    # with keys the basename and values the normalised array
    if isinstance(array, dict):
        colors = get_N_HexCol(len(array.keys()))
        # Remove _1 or _2 etc to get cell_condition-protein
        print 'Next come array.keys'
        print sorted(array.keys())
        proteins = [ re.sub(r'_\d{1}.*\\.$', '\\.', k.split('.')[0]) for k in sorted(array.keys()) ]
        print 'Next come proteins'
        print proteins

        # create a dictionary to hold colors used in the average plot
        # All replicates have the same color but different shade
        colors_per_protein = {}
        counter = 0
        for p in set(proteins):
            reps = [x for x, library in enumerate(proteins) if p in library]

            if len(reps) == 1:
                colors_per_protein[sorted(array.keys())[x]] = colors[counter]
            else:
                offset = 1
                for x in reps:
                    colors_per_protein[sorted(array.keys())[x]] = color_variant(colors[counter], offset)
                    # offset smaller than 50 doesn't allow distinction between replicates
                    offset += 50

            counter += 1

        for key in sorted(array.keys()):

            metaseq.plotutils.ci_plot(
                    window,
                    array.get(key),
                    ax=ax,
                    line_kwargs = dict(color=colors_per_protein.get(key), label = key),
                    fill_kwargs = dict(color=colors_per_protein.get(key), alpha=0.3)
                    )
    elif name is not None:
        ax.plot(window, array.mean(axis=0), color = 'r', label = name)
        # Add a vertical line at TSS/gene-start
        ax.axvline(0, color='k', linestyle='--')

    else:
        print "You did not pass a name for the sample."
    
    # Add labels and legend
    if feature_type == 'genes':
        ax.set_xlabel('Percentage of ' + feature_type)
    else:
        ax.set_xlabel('Distance from ' + feature_type +' (bp)')
    ax.set_ylabel('Normalised average read coverage (per M mapped reads)')
    ax.legend(loc = 2, frameon = False, fontsize = 14, labelspacing = .3, handletextpad = 0.2)
    
    pp.savefig(fig)
    plt.close(fig)
    # }}}

def plot_tss(array, window, name, pp, subsets = None, features = None, feature_type = "Genes", expression = None):# {{{

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
            'sort_by':array.mean(axis=1)
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
    elif expression is None:
        arguments['sort_by'] = array.mean(axis=1)
        print "No subsets provided"# }}}
    
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

    if feature_type == 'genes':
        arguments['vmin'] = 1
        arguments['vmax'] = 99
        arguments['imshow_kwargs'] = dict(interpolation = 'none')

    fig = metaseq.plotutils.imshow(array, **arguments)

    if features is not None and expression is not None:
        fig = format_axes(fig, feature_type = feature_type, subsets = True, by = arguments['subset_by'], order = arguments['subset_order'])
    elif features is not None and subsets is not None:
        fig = format_axes(fig, feature_type = feature_type, subsets = True, by = arguments['subset_by'], order = arguments['subset_order'])
    else:
        fig = format_axes(fig, feature_type = feature_type)

    pp.savefig(fig)
    plt.close(fig)
    # }}}

# }}}

# # For each bam calculate signal and plot it# {{{
# # The windows we'll get signal over
def main():
    from matplotlib.backends.backend_pdf import PdfPages
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
    per_feature_type = {'genes':dict(), 'tss':dict(), 'gene_start':dict()}
    
    # One file for all replicates in a given chip# {{{
    for x in names:
        # Create pdfpage object
        pp = PdfPages(plot_path + x + '-profiles.pdf')

        npz = [ k for k in normalised_arrays.keys() if re.match(x,k) ]
        
        for n in npz:
            print 'Plotting for %s' % n
            # This should only return one file
            features = [ os.path.join(args.path, f) for f in os.listdir(args.path) if re.match(re.escape(n) + r'.features',f) ][0]
            features = pybedtools.BedTool(features)

            window = os.path.splitext(n)[1].replace(".", "")
            pattern = re.compile('(\d+)-(\w+)-(\d+)')
            # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
            matching = pattern.match(window)
            
            if matching:
                upstream, feature_type, downstream = matching.groups()
                # bins are of size 100bp in count-tags.py
                x = np.linspace(-int(upstream), int(downstream), 100)
                
                if feature_type == 'tss':
                    per_feature_type['tss'][n] = normalised_arrays.get(n)
                elif feature_type == 'gene_start':
                    per_feature_type['gene_start'][n] = normalised_arrays.get(n)
            elif window == 'genes':
                feature_type = 'genes'
                # gene array goes from o) to 100% in bins of 100bp
                x = np.linspace(0, 100, 100)

                per_feature_type['genes'][n] = normalised_arrays.get(n)

            else:
                print 'Unknown feature type for %s. Exiting now.' % n
                pp.close()
                sys.exit()
            
            if feature_type == 'gene_start' or feature_type == 'genes':
                plot_tss(array = normalised_arrays.get(n), window = x, name=n,
                        pp = pp, features = features, feature_type = feature_type)
                if args.subsets:
                    plot_tss(array = normalised_arrays.get(n), window = x, name=n,
                            pp = pp, features = features, feature_type = feature_type,
                            subsets = args.subsets)
                    if args.expression:
                        plot_tss(array = normalised_arrays.get(n), window = x, name=n,
                                pp = pp, features = features, feature_type = feature_type,
                                expression = args.expression)
        pp.close()
# }}}

    # Plotting all ChIP experiments together.# {{{
    # TODO: between sample normalisation?
    if args.matching == "*":
        identifier = strftime("%Y-%m-%d-%Hh%Mm", gmtime())
    else:
        identifier = args.matching.replace('*', '')

    pp = PdfPages(plot_path + identifier + '.all-profiles.pdf')
    print 'Printing all profiles'

    for key in per_feature_type.keys():
        # key %in% tss, gene_start, genes
        print key
        window = [ w.split('.')[1] for w in per_feature_type.get(key).keys() ]
        print per_feature_type.get(key).keys()
        if len(set(window)) == 1:
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
            plot_average(per_feature_type[key], x, pp = pp, feature_type = key)
        else:
            print 'Can not plot experiment together. Tags calculated on different windows'
            
    pp.close()
# }}}

# }}}

main()
