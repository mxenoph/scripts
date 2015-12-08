#!/usr/bin/env python

# ---------------------------------------------------------
# Calculate tags around TSS/gene/region and normalise to 
# input and library size.
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
# Required arguments
parser.add_argument('-i', '--ip', nargs = '+', metavar = "file", type = str, required = True, help = 'One or more treated(IP) bam files')
parser.add_argument('-c', '--ctrl', nargs = '+', metavar = "file", type = str, required = True, help = 'Control bam files')
parser.add_argument('-o', '--out_dir', metavar ="path", type = str, required = True)
parser.add_argument('--fs', required = False, default = 200, help = "Fragment size from library preparation")
parser.add_argument('-u', '--upstream', required = False, default = 2000, help = "bp upstream of tss/peak centre")
parser.add_argument('-d', '--downstream', required = False, default = 500, help ="bp downstream of tss/peak centre")
parser.add_argument('-g', '--gtf', type = str,
        default = '/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.70.gtf',
        help= 'GTF for the ensembl annotation.')
parser.add_argument('-r', '--regions', type = str, help= 'BED file containing regions of interest.')
parser.add_argument('-f', '--full', default = True, help="If set to true then count tags for gene/peak")

args = parser.parse_args()

if not args.out_dir.endswith(os.sep):
    args.out_dir = args.out_dir + os.sep

args.out_dir = args.out_dir + 'profiles/'

# Create plot directory if it doesn't exist
if not os.path.exists(args.out_dir):
    os.makedirs(args.out_dir)

pattern = re.compile(".*[I|i]nput.*")
# Remove any control files passed as ip accidentally
args.ip = [x for x in args.ip if not pattern.match(x)]
# Remove any ip files passed as control accidentally
#args.ctrl = [x for x in args.ctrl if pattern.match(x)]

# Create empty array/list
basename = []
for i in args.ip:
    base = os.path.splitext(os.path.basename(i))[0]
    basename.append(base)
#}}}

# Functions# {{{
class ref:
    def __init__(self, obj): self.obj = obj
    def get(self):    return self.obj
    def set(self, obj):      self.obj = obj

def genes_generator(database):# {{{
    """
    Generator function to yield TSS +/- Kb of each annotated gene
    """
    for gene in database.features_of_type('gene'):
        if (re.match('chr', gene.chrom)):
            yield asinterval(gene)
# }}}

def gene_start_generator(database):# {{{
    """
    Generator function to yield TSS +/- Kb of each annotated gene
    """
    for gene in database.features_of_type('gene'):
        if (re.match('chr', gene.chrom)):
            yield TSS(asinterval(gene), upstream = args.upstream, downstream = args.downstream)
# }}}

def tss_generator(database):# {{{
    """
    Generator function to yield TSS +/- kb of each annotated transcript
    """
    for transcript in database.features_of_type('transcript'):
        if (re.match('chr', transcript.chrom)):
            yield TSS(asinterval(transcript), upstream = args.upstream, downstream = args.downstream)# }}}

# Create arrays in parallel, and save to disk for later # {{{
def count_tags (ip, ctrl, features, description):
    "This counts mapped reads for ip and input and normalizes them by library size and million mapped reads"
    from metaseq import persistence
    import multiprocessing
    processes = multiprocessing.cpu_count()

    # Use prefix of ip file as filename for numpy array (npz)
    basename = args.out_dir + os.path.splitext(os.path.basename(ip))[0]
    # description is a description of the features provided
    basename = basename + '.' + description

    output = basename + '.npz'
    print output
    sys.exit
    # Run if file does not exist and experiment has no replicates
    if os.path.exists(output):
        print ".npz already exists; using that"
    else:
        # Read in bam files
        ip = metaseq.genomic_signal(ip, 'bam')
        ctrl = metaseq.genomic_signal(ctrl, 'bam')
        
        # Create arrays in parallel
        print "Building the IP array for %s using %s processors" % (basename, processes)
        #ip_array = ip.array(features, bins = 100, processes = processes, fragment_size = int(args.fs))
        ip_array = ip.array(features, bins = 100, processes = processes, shift_width = int(args.fs)/2)
        print "Building the input array for %s using %s processors" % (basename, processes)
        #ctrl_array = ctrl.array(features, bins = 100, processes = processes, fragment_size = int(args.fs))
        ctrl_array = ctrl.array(features, bins = 100, processes = processes, shift_width = int(args.fs)/2)
        
        # Normalize to library size
        ip_array /= ip.mapped_read_count() / 1e6
        ctrl_array /= ctrl.mapped_read_count() / 1e6

        # Cache to disk (will be saved as "example.npz" and "example.features")
        persistence.save_features_and_arrays(
                features = features,
                arrays={'ip': ip_array, 'input': ctrl_array},
                prefix = basename,
                link_features = True,
                overwrite = True)
        print "done"
    return;# }}}

# # Create database from ensembl GTF if it does not already exist# {{{
def create_features():
    gff_filename = args.gtf
    db_filename = gff_filename + '.db'
    
    if not os.path.exists(db_filename):
        gffutils.create_db(gff_filename, db_filename)
        
    db = gffutils.FeatureDB(db_filename)
    
    suffix = '.metaseq.' + str(args.upstream) + '-tss-' + str(args.downstream) + '.gtf'
    filename = gff_filename.replace('.gtf', suffix)
    # Here we only create if needed, caching to disk.
    if not os.path.exists(filename):
        print "Creating tss: %s " % filename
        # A BedTool made out of a generator, and saved to file.
        tss = pybedtools.BedTool(tss_generator(db)).saveas(filename)
    else:
        tss = pybedtools.BedTool(filename)
        
    filename = gff_filename.replace('.gtf', '.metaseq.genes.gtf')
    if not os.path.exists(filename):
        # A BedTool made out of a generator, and saved to file.
        genes = pybedtools.BedTool(gff_filename)\
                        .saveas(filename)
    else:
        genes = pybedtools.BedTool(filename)
                
    suffix = '.metaseq.' + str(args.upstream) + '-gene_start-' + str(args.downstream) + '.gtf'
    filename = gff_filename.replace('.gtf', suffix)
    if not os.path.exists(filename):
        # A BedTool made out of a generator, and saved to file.
        gene_start = pybedtools.BedTool(gene_start_generator(db))\
                .saveas(filename)
    else:
        gene_start = pybedtools.BedTool(filename)
                
    #tsses_1kb = tsses.slop(b=1000, genome='mm10', output = gff_filename.replace('.gtf', '.tss_1kb.gtf'))
    return {'tss':tss, 'gene_start':gene_start, 'genes':genes}
# }}}

# }}}

# For each bam calculate signal and plot it# {{{
def main():
    # Create genomic_signal objects that point to data files
    for i in range(len(args.ip)):
        if len(args.ip) == len(args.ctrl):
            bams = {'ip':args.ip[i], 'ctrl':args.ctrl[i]}
            print 'Same number of IP and input BAM bams provided. Assuming consistency in listing bams.'
            print 'Counting tags for IP: %(ip)s and Input: %(ctrl)s' % bams
        elif len(args.ctrl) == 1:
            bams = {'ip':args.ip[i], 'ctrl':args.ctrl[0]}
            print 'Only one input BAM file provided; this is used to correct background noise for all IP BAM bams.'
            print 'Counting tags for IP: %(ip)s and Input: %(ctrl)s' % bams
        else:
            print 'Different number of IP and input BAM bams provided. I do not know how these experiments are associated.'

        features = create_features()
        count_tags(ip = bams['ip'], ctrl = bams['ctrl'],\
                features = features['gene_start'], description = str(args.upstream) + '-gene_start-' + str(args.downstream))

        count_tags(ip = bams['ip'], ctrl = bams['ctrl'],\
                features = features['tss'], description = str(args.upstream) + '-tss-' + str(args.downstream))

        count_tags(ip = bams['ip'], ctrl = bams['ctrl'],\
                features = features['genes'], description = 'genes' )
        
        if not args.regions:
            print 'Regions not provided'
        else:
            'ToDo: read in the peaks and do similar analysis to genes'
# }}}

main()
