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
import pybedtools
import metaseq
from pylab import *
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval# }}}

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--ip', nargs='+', metavar="file", type=str, required=True, help= 'One or more treated(IP) bam files')
parser.add_argument('-c', '--ctrl', nargs='+', metavar="file", type=str, required=True, help= 'Control bam files')
parser.add_argument('-a', '--assembly', type=str, default='mm9', help= 'Assembly for the ensembl annotation. Default = mm9')
parser.add_argument('-f', '--gtf', metavar="file", type=str, required=False, help= 'Alternative feature gtf')
parser.add_argument('-e', '--expr', metavar="file", type=str, required=False, help= 'Deseq output file')
parser.add_argument('-b', '--bound', metavar="file", type=str, required=False, help= 'TSV file with gene_id (ensembl) and bound (1/0) fields')
parser.add_argument('-l', '--label', type=str, required=False, default= strftime("%Y-%m-%d-%Hh%Mm", gmtime()), help= 'Label to be used in plottting output file')
parser.add_argument('-o', '--out_dir', metavar="path", type=str, required=True)

args = parser.parse_args()

if not args.out_dir.endswith(os.sep):
    args.out_dir = args.out_dir + os.sep

plot_dir = args.out_dir + 'plots/'
# Create plot directory if it doesn't exist
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

pattern = re.compile(".*[I|i]nput.*")
# Remove any control files passed as ip accidentally
args.ip = [x for x in args.ip if not pattern.match(x)]
# Remove any ip files passed as control accidentally
args.ctrl = [x for x in args.ctrl if pattern.match(x)]

# Create empty array/list
basename = []
for i in args.ip:
    base = os.path.splitext(os.path.basename(i))[0]
    basename.append(base)#}}}

# Functions# {{{
class ref:
    def __init__(self, obj): self.obj = obj
    def get(self):    return self.obj
    def set(self, obj):      self.obj = obj

def tss_generator():# {{{
    """
    Generator function to yield TSS +/- 1kb of each annotated gene
    """
    for transcript in db.features_of_type('gene'):
        if (re.match('chr', transcript.chrom)):
            yield TSS(asinterval(transcript), upstream=1000, downstream=1000)# }}}

# Create arrays in parallel, and save to disk for later
def calc_signal ( ip, ctrl, anchor, basename):# {{{
    "This counts mapped reads for ip and input and normalizes them by library size and million mapped reads"
    from metaseq import persistence
    import multiprocessing
    processes = multiprocessing.cpu_count()
    
    out = basename + '.npz'
    # Run if file does not exist and experiment has no replicates
    if not os.path.exists(out):
        # Create arrays in parallel
        ip_array = ip_signal.array(anchor, bins=100, processes=processes)
        input_array = input_signal.array(anchor, bins=100, processes=processes)
        
        # Normalize to library size
        ip_array /= ip_signal.mapped_read_count() / 1e6
        input_array /= input_signal.mapped_read_count() / 1e6
        
        # Cache to disk (will be saved as "example.npz" and "example.features")
        persistence.save_features_and_arrays(
                features=anchor,
                arrays={'ip': ip_array, 'input': input_array},
                prefix=basename,
                link_features=True,
                overwrite=True)
    return;# }}}

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

def plot_norm_signals (norm_sub, label, window):# {{{
    # Create a meaningful x-axis
    import numpy as np
    #x = np.linspace(-1000, 1000, 100)
    window = window / 2
    x = np.linspace(-window, window, 100)
    
    # Initial plot of average signal over TSSs
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    n_samples = len(norm_sub.keys())
    colors = get_N_HexCol(n_samples)
    
    for key in norm_sub.keys():
        ax.plot(x, norm_sub.get(key).mean(axis=0), color=colors[norm_sub.keys().index(key)], label=key)
    
    # Add a vertical line at the TSS
    ax.axvline(0, linestyle=':', color='k')
    ax.axvline(75, linestyle=':', color='c')
    ax.axvline(150, linestyle=':', color='r')
    ax.axvline(-75, linestyle=':', color='c')
    ax.axvline(-150, linestyle=':', color='r')
    
    # Add labels and legend
    ax.set_xlabel('Distance from ' + label +' (bp)')
    ax.set_ylabel('Normalised average read coverage (per million mapped reads)')
    ax.legend(loc=2, frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)

    return fig;# }}}

def signal_expr(norm_sub, expr, window, label = ' (DE)'):# {{{
    import numpy as np
    from matplotlib import pyplot as plt

    window = window / 2
    x = np.linspace(-window, window, 100)
    
    de = plt.figure()
    # Depends on the replicate number
    number_of_subplots = len(norm_sub.keys())

    def plot_by_expression(normalized_subtracted, expr, axes, title):# {{{
        metaseq.plotutils.ci_plot(
                x,
                normalized_subtracted[((expr.log2FoldChange > 0) & (expr.padj <= 0.05)).values, :],
                line_kwargs=dict(color='#fe9829', label='up'),
                fill_kwargs=dict(color='#fe9829', alpha=0.3),
                ax= axes)
        metaseq.plotutils.ci_plot(
                x,
                normalized_subtracted[((expr.log2FoldChange < 0) & (expr.padj <= 0.05)).values, :],
                line_kwargs=dict(color='#8e3104', label='down'),
                fill_kwargs=dict(color='#8e3104', alpha=0.3),
                ax= axes)
        metaseq.plotutils.ci_plot(
                x,
                normalized_subtracted[((expr.padj > 0.05)).values, :],
                line_kwargs=dict(color='.5', label='unchanged'),
                fill_kwargs=dict(color='.5', alpha=0.3),
                ax= axes)
        
        axes.set_ylabel('Average\nenrichment')
        axes.set_title(title)
        axes.set_xlim([-window, window])
        axes.axvline(0, linestyle=':', color='k')
        axes.legend(loc=2, frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
        
        return axes;# }}}

    for key in norm_sub.keys():
        normalized_subtracted = norm_sub.get(key)
        # This will return the index number for the key which will then be used to define subplot
        v = norm_sub.keys().index(key) + 1
        import math
        v = int(str(int(math.ceil(number_of_subplots/2))) + str(2) + str(v))
      
        ax = de.add_subplot(v)
        title = key + label
        #axes = plot_by_expression(normalized_subtracted, expr, ax, title)
        plot_by_expression(normalized_subtracted, expr, ax, title)
        
    de.tight_layout()
    return de;
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
    x = np.linspace(-window, window, 100)
    
    n_samples = len(norm_sub.keys())
    colors = get_N_HexCol(n_samples)

    fig = plt.figure()
    # For all bound TSS
    ax = fig.add_subplot(111)
    heatmaps = dict()
    norm_sub_bound = dict()
    
    # Converting numbers to strings and making a matrix to hold the information for bound/not# {{{
    subset_by = [str(num) for num in subset.boolean.values]
    subset_by = [w.replace('1', 'Bound') for w in subset_by]
    subset_by = [w.replace('0', 'Not bound') for w in subset_by]
    subset_by = array(subset_by, dtype="|S32")

    subset_order = ['Not bound', 'Bound']
    # }}}

    for key in norm_sub.keys():
        normalized_subtracted = norm_sub.get(key)

        heatmaps[key]= metaseq.plotutils.imshow(normalized_subtracted,
                x,
                features=tsses,
                vmin=5, vmax=95, percentile=True,
                sort_by=normalized_subtracted.mean(axis=1),
                subset_by= subset_by,
                subset_order=subset_order,
                line_kwargs=[dict(color='k', label='Not bound'), dict(color='r', label='Bound')],
                fill_kwargs=[dict(color='k', alpha=0.3), dict(color='r', alpha=0.3)])
        # adding labels to heatmap
        metaseq.plotutils.add_labels_to_subsets(heatmaps[key].array_axes, subset_by=subset_by, subset_order=subset_order,)
        # adding labels to average signal plot
        heatmaps[key].line_axes.legend(loc='best', frameon=False)
        heatmaps[key].suptitle(key + "\n(Bound genes: " + str(len(subset[subset_by=='Bound'])) + ' )')
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

# }}}

# Create database from ensembl GTF if it does not already exist# {{{
import pandas as df
annotation_info = df.read_csv('/nfs/research2/bertone/user/mxenoph/genome_dir/assemblies-annotations.config', sep="\t")
# Otherwise pandas subsetting returns a series (despite only one string there) instead of a single string
gff_filename = annotation_info.loc[annotation_info.assembly == args.assembly].annotation.iloc[0]

#gff_filename = '/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/Mus_musculus.NCBIM37.67_conv.gtf'
db_filename = gff_filename + '.db'

if not os.path.exists(db_filename):
    gffutils.create_db(gff_filename, db_filename)

db = gffutils.FeatureDB(db_filename)# }}}

# Create a GTF for TSS from gene start for all genes on chromosomes# {{{
tss_filename= gff_filename.replace('.gtf', '.tss.gtf')
# At the moment have created the tss gtf by grep chr on the gtf creating from the script but must be a way to do
# it with gffutils
tss_filename= gff_filename.replace('.gtf', '.tss_nocontig.gtf')

# Here we only create if needed, caching to disk.
if not os.path.exists(tss_filename):
    # A BedTool made out of a generator, and saved to file.
    tsses = pybedtools.BedTool(tss_generator()).saveas(tss_filename)
    # }}}

# For each bam calculate signal and plot it# {{{
# The windows we'll get signal over
tsses = pybedtools.BedTool(tss_filename)
# for looking at other features

# If alternative features passed then create the object
if args.gtf:
    alternative_features = pybedtools.BedTool(args.gtf)

norm_sub = dict()
norm_sub_alt_feat = dict()

# Create genomic_signal objects that point to data files
for i in range(len(args.ip)):
    base = os.path.splitext(os.path.basename(args.ip[i]))[0]
    ip_filename = args.ip[i]
    input_filename = args.ctrl[i]
    print ip_filename
    print input_filename

    ip_signal = metaseq.genomic_signal(ip_filename, 'bam')
    input_signal = metaseq.genomic_signal(input_filename, 'bam')

    calc_signal(ip_signal, input_signal, tsses, base)

    # Load the windows and arrays
    from metaseq import persistence
    features, arrays = persistence.load_features_and_arrays(prefix=base)
    
    # Normalize IP to the control
    normalized_subtracted = arrays['ip'] - arrays['input']
    norm_sub[ip_filename] = normalized_subtracted

    try:# {{{
        alternative_features
    except:
        print "Only computing enrichment around the TSS."
    else:
        print 'Computing for alternative features'
        gtf_base = base + '.alternative_features'

        from joblib import Parallel, delayed
        import multiprocessing
        inputs = range(1, len(alternative_features))
        processes = multiprocessing.cpu_count()
        
        # returns the width of a range (i) in the pybedtool object
        def interval_width(i, bedtool):
            return len(bedtool[i])

        widths = Parallel(n_jobs = processes)(delayed(interval_width)(i, alternative_features) for i in inputs)
        window = most_common(widths)
        
        calc_signal(ip_signal, input_signal, alternative_features, gtf_base)
        features, arrays = persistence.load_features_and_arrays(prefix = gtf_base)
        # Normalize IP to the control
        normalized_subtracted = arrays['ip'] - arrays['input']
        norm_sub_alt_feat[ip_filename] = normalized_subtracted # }}}

# }}}

# Plotting ## {{{
from matplotlib.backends.backend_pdf import PdfPages
# Create pdfpage object
pp = PdfPages(plot_dir + args.label + '-averageSignal.pdf')
pp.savefig(plot_norm_signals(norm_sub, 'TSS', 2000))

try:# {{{
    alternative_features
except:
    print "Only computing enrichment around the TSS."
else:
    from joblib import Parallel, delayed
    import multiprocessing
    inputs = range(1, len(alternative_features))
    processes = multiprocessing.cpu_count()
    
    # returns the width of a range (i) in the pybedtool object
    def interval_width(i, bedtool):
        return len(bedtool[i])

    widths = Parallel(n_jobs = processes)(delayed(interval_width)(i, alternative_features) for i in inputs)
    window = most_common(widths)
    
    pp.savefig(plot_norm_signals(norm_sub_alt_feat, 'midpoint', window)) # }}}
    
if args.expr:
    from metaseq.results_table import DESeqResults
    expr = DESeqResults(args.expr, import_kwargs=dict(index_col=0))
    expr = expr.reindex_to(tsses, attribute='gene_id')
    pp.savefig(signal_expr(norm_sub, expr, 2000))

if args.bound:
    from metaseq.results_table import ResultsTable
    bound = ResultsTable(args.bound, import_kwargs=dict(index_col=0))
    # reindexing the bound dataframe 
    bound = bound.reindex_to(tsses, attribute = 'gene_id')
    # Pass pp object to function and print within function
    signal_bound(norm_sub, bound, 2000, pp, expr=expr)

pp.close()
# }}}
