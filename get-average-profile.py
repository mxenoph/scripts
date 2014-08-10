#!/usr/bin/python

# ---------------------------------------------------------
# Get TSS +/- 1kb for all annotated transcripts.
#

# Import modules# {{{
import os, sys, argparse, re
import gffutils
import pybedtools
import metaseq
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval# }}}

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--ip', nargs='+', metavar="file", type=str, required=True, help= 'One or more treated(IP) bam files')
parser.add_argument('-c', '--ctrl', metavar="file", type=str, required=True, help= 'Control bam file')
parser.add_argument('-e', '--expr', metavar="file", type=str, required=False, help= 'Deseq output file')
parser.add_argument('-o', '--out_dir', metavar="path", type=str, required=True)

args = parser.parse_args()

if not args.out_dir.endswith(os.sep):
    args.out_dir = args.out_dir + os.sep

plot_dir = args.out_dir + 'plots/'
# Create plot directory if it doesn't exist
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

pattern = re.compile(".*[I|i]nput.*")
args.ip = [x for x in args.ip if not pattern.match(x)]

# Create empty array/list
basename = []
for i in args.ip:
    base = os.path.splitext(os.path.basename(i))[0]
    basename.append(base)#}}}

# Functions# {{{
# Create arrays in parallel, and save to disk for later
def calc_signal ( ip, ctrl, anchor, basename ):
    "This counts mapped reads for ip and input and normalizes them by library size and million mapped reads"
    from metaseq import persistence
    import multiprocessing
    processes = multiprocessing.cpu_count()
    
    out = basename + '.npz'
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
    return;

def plot_signals (arrays, base):
    # Create a meaningful x-axis
    import numpy as np
    x = np.linspace(-1000, 1000, 100)
    
    # Initial plot of average signal over TSSs
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, arrays['ip'].mean(axis=0), color='r', label=base)
    ax.plot(x, arrays['input'].mean(axis=0), color='k', label='input')
    
    # Add a vertical line at the TSS
    ax.axvline(0, linestyle=':', color='k')
    
    # Add labels and legend
    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best');

    return fig;

def heatmap_signal_norm (normalized_subtracted):
    # Create a meaningful x-axis
    import numpy as np
    x = np.linspace(-1000, 1000, 100)
    ## First version of a plot that includes a heatmap of the array
    from matplotlib import pyplot as plt
    plt.rcParams['font.size'] = 8
    
    fig = metaseq.plotutils.imshow(
            normalized_subtracted, # The array to plot
            x=x, # X-axis to use
            # Make the colorbar limits go from 5th to 99th percentile. 
            # `percentile=True` means treat vmin/vmax as percentiles rather than
            # actual values.
            vmin=5, vmax=99,  percentile=True,
            line_kwargs=dict(color='k', label='All'), # Style for the average line plot
            fill_kwargs=dict(color='k', alpha=0.3), # Style for the +/- 95% CI band surrounding the average line
            sort_by=normalized_subtracted.mean(axis=1) # Additionally, sort by mean signal
            )
    
    # Label axes, add dotted lines indicating TSS
    fig.line_axes.set_ylabel('Average enrichment');
    fig.line_axes.set_xlabel('Distance from TSS (bp)');
    
    fig.array_axes.set_ylabel('Genes')
    fig.array_axes.set_xticklabels([])
    
    fig.array_axes.axvline(0, linestyle=':', color='k')
    fig.line_axes.axvline(0, linestyle=':', color='k')

    return fig;

#def heatmap_signal_expr (normalized_subtracted, expr):
def heatmap_signal_expr (normalized_subtracted, expr, fc_order):
    # Create a meaningful x-axis
    import numpy as np
    x = np.linspace(-1000, 1000, 100)

    from matplotlib import pyplot as plt
    fig = metaseq.plotutils.imshow(
            normalized_subtracted[fc_order],
            x=x,
            vmin=5, vmax=99,  percentile=True,
            line_kwargs=dict(color='k', label='All'),
            fill_kwargs=dict(color='k', alpha=0.3),
            #sort_by=normalized_subtracted.mean(axis=1),
            # Additionally specify height_ratios:
            height_ratios=(3, 1, 1)
            )
    
    # `fig.gs` contains the `matplotlib.gridspec.GridSpec` object,
    # so we can now create the new axes.
    bottom_axes = plt.subplot(fig.gs[2, 0])

    # Signal over TSSs of transcripts that were activated upon knockdown.
    metaseq.plotutils.ci_plot(
            x,
            normalized_subtracted[((expr.log2FoldChange > 0) & (expr.padj <= 0.01)).values, :],
            line_kwargs=dict(color='#fe9829', label='up'),
            fill_kwargs=dict(color='#fe9829', alpha=0.3),
            ax=bottom_axes)
    # Signal over TSSs of transcripts that were repressed upon knockdown
    metaseq.plotutils.ci_plot(
            x,
            normalized_subtracted[((expr.log2FoldChange < -1) & (expr.padj <= 0.01)).values, :],
            line_kwargs=dict(color='#8e3104', label='down'),
            fill_kwargs=dict(color='#8e3104', alpha=0.3),
            ax=bottom_axes)

    # Signal over TSSs tof transcripts that did not change upon knockdown
    metaseq.plotutils.ci_plot(
            x,
            normalized_subtracted[((expr.padj > 0.01)).values, :],
            line_kwargs=dict(color='.5', label='unchanged'),
            fill_kwargs=dict(color='.5', alpha=0.3),
            ax=bottom_axes)
    
    # Clean up redundant x tick labels, and add axes labels
    fig.line_axes.set_xticklabels([])
    fig.array_axes.set_xticklabels([])
    fig.line_axes.set_ylabel('Average\nenrichement')
    fig.array_axes.set_ylabel('Genes')
    bottom_axes.set_ylabel('Average\nenrichment')
    bottom_axes.set_xlabel('Distance from TSS (bp)')
    fig.cax.set_ylabel('Enrichment')
    
    # Add the vertical lines for TSS position to all axes
    for ax in [fig.line_axes, fig.array_axes, bottom_axes]:
        ax.axvline(0, linestyle=':', color='k')

        # Nice legend
        bottom_axes.legend(loc='best', frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)
        fig.subplots_adjust(left=0.3, right=0.8, bottom=0.05)

    return fig;

# }}}

# Create database from ensembl GTF if it does not already exist# {{{
gff_filename = '/nfs/research2/bertone/user/mxenoph/genome_dir/M_musculus_9/MM9.maps/Mus_musculus.NCBIM37.67_conv.gtf'
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
    def tss_generator():
        """
        Generator function to yield TSS +/- 1kb of each annotated gene
        """
        for transcript in db.features_of_type('gene'):
            if (re.match('chr', transcript.chrom)):
                yield TSS(asinterval(transcript), upstream=1000, downstream=1000)
                
    # A BedTool made out of a generator, and saved to file.
    tsses = pybedtools.BedTool(tss_generator()).saveas(tss_filename)
    # }}}

# For each bam calculate signal and plot it# {{{
# The windows we'll get signal over
tsses = pybedtools.BedTool(tss_filename)


# Create genomic_signal objects that point to data files
for i in args.ip:
    base = os.path.splitext(os.path.basename(i))[0]
    ip_filename = i
    input_filename = args.ctrl

    #ip_filename= '/nfs/research2/bertone/user/mxenoph/hendrich/chip/factor_2014/mm9/bowtie/mEpiSC_H3K4me3.bam'
    #input_filename='/nfs/research2/bertone/user/mxenoph/hendrich/chip/factor_2014/mm9/bowtie/mEpiSC_input.bam'
    ip_signal = metaseq.genomic_signal(ip_filename, 'bam')
    input_signal = metaseq.genomic_signal(input_filename, 'bam')

    calc_signal(ip_signal, input_signal, tsses, base)
    
    # Load the windows and arrays
    from metaseq import persistence
    features, arrays = persistence.load_features_and_arrays(prefix=base)
    
    # Normalize IP to the control
    normalized_subtracted = arrays['ip'] - arrays['input']
    
    from matplotlib.backends.backend_pdf import PdfPages
    # Create pdfpage object
    pp = PdfPages(plot_dir + base + '-averageSignal.pdf')

    pp.savefig(plot_signals(arrays, base))
    pp.savefig(heatmap_signal_norm(normalized_subtracted))

    if args.expr: #'expr' in args:
        from metaseq.results_table import ResultsTable
        from metaseq.results_table import DESeqResults
        expr = ResultsTable(args.expr, import_kwargs=dict(index_col=0))
        expr = expr.reindex_to(tsses, attribute='gene_id')
        test= DESeqResults(args.expr, import_kwargs=dict(index_col=0))
        test = test.reindex_to(tsses, attribute='gene_id')
        up = test.upregulated(thresh=0.05, idx=True, col='padj').ravel().nonzero()[0]
        down = test.downregulated(thresh=0.05, idx=True, col='padj').ravel().nonzero()[0]
        unch = test.unchanged(thresh=0.05, idx=True, col='padj').ravel().nonzero()[0]
        import numpy as np
        fc_order = np.concatenate([up, unch , down])
        pp.savefig(heatmap_signal_expr(normalized_subtracted, expr, fc_order))

pp.close()# }}}
