#!/usr/bin/env python

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
parser.add_argument('-c', '--ctrl', nargs='+', metavar="file", type=str, required=True, help= 'Control bam files')
parser.add_argument('-a', '--assembly', type=str, default='mm9', help= 'Assembly for the ensembl annotation. Default = mm9')
parser.add_argument('-f', '--gtf', metavar="file", type=str, required=False, help= 'Alternative feature gtf')
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

def tss_generator():# {{{
    """
    Generator function to yield TSS +/- 1kb of each annotated gene
    """
    for transcript in db.features_of_type('gene'):
        if (re.match('chr', transcript.chrom)):
            yield TSS(asinterval(transcript), upstream=1000, downstream=1000)# }}}

# Create arrays in parallel, and save to disk for later
def calc_signal ( ip, ctrl, anchor, basename, replicates =False ):# {{{
    "This counts mapped reads for ip and input and normalizes them by library size and million mapped reads"
    from metaseq import persistence
    import multiprocessing
    processes = multiprocessing.cpu_count()
    
    out = basename + '.npz'
    # Run if file does not exist and experiment has no replicates
    if not os.path.exists(out) & (not replicates):
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

def plot_norm_signals (norm_sub):# {{{
    # Create a meaningful x-axis
    import numpy as np
    x = np.linspace(-1000, 1000, 100)
    
    # Initial plot of average signal over TSSs
    from matplotlib import pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    colors = ['b', 'g', 'r', 'c', 'k', 'm']
    colors = get_N_HexCol(12)
    
    for key in norm_sub.keys():
        ax.plot(x, norm_sub.get(key).mean(axis=0), color=colors[norm_sub.keys().index(key)], label=key)
    
    # Add a vertical line at the TSS
    ax.axvline(0, linestyle=':', color='k')
    ax.axvline(75, linestyle=':', color='c')
    ax.axvline(150, linestyle=':', color='r')
    ax.axvline(-75, linestyle=':', color='c')
    ax.axvline(-150, linestyle=':', color='r')
    
    # Add labels and legend
    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Normalised average read coverage (per million mapped reads)')
    ax.legend(loc=2, frameon=False, fontsize=8, labelspacing=.3, handletextpad=0.2)

    return fig;# }}}
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

    try:
        alternative_features
    except:
        print "Only computing enrichment around the TSS."
    else:
        calc_signal(ip_signal, input_signal, alternative_features, base + '.mb3_peak')


    # Load the windows and arrays
    from metaseq import persistence
    features, arrays = persistence.load_features_and_arrays(prefix=base)
    
    # Normalize IP to the control
    normalized_subtracted = arrays['ip'] - arrays['input']
    norm_sub[ip_filename] = normalized_subtracted


from matplotlib.backends.backend_pdf import PdfPages
# Create pdfpage object
pp = PdfPages(plot_dir + 'm2-mi2b' + '-averageSignal.pdf')
pp.savefig(plot_norm_signals(norm_sub))
pp.close()
# }}}
