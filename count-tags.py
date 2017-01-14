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
from pybedtools.featurefuncs import less_than
from pybedtools.featurefuncs import greater_than
from gffutils.helpers import asinterval# }}}

# Parsing command line arguments and creating output subdirectories# {{{
parser = argparse.ArgumentParser()
# Required arguments
parser.add_argument('-i', '--ip', nargs = '+', metavar = "file", type = str, required = True, help = 'One or more treated(IP) bam files, or bigwig files if -b flag is set')
parser.add_argument('-c', '--ctrl', nargs = '+', metavar = "file", type = str, required = False,  help = 'Control bam files')
parser.add_argument('-o', '--out_dir', metavar ="path", type = str, required = True)
parser.add_argument('--fs', required = False, default = 200, help = "Fragment size from library preparation")
parser.add_argument('-u', '--upstream', required = False, default = 2000, help = "bp upstream of tss/peak centre")
parser.add_argument('-d', '--downstream', required = False, default = 500, help ="bp downstream of tss/peak centre")
parser.add_argument('-g', '--gtf', type = str,
        default = '/nfs/research2/bertone/user/mxenoph/common/genome/MM10/Mus_musculus.GRCm38.75.gtf',
        help= 'GTF for the ensembl annotation.')
parser.add_argument('-a', '--assembly', type = str, default = 'mm10', help= 'Assembly (default:mm10). Used to retrieve chromosomes lengths from ucsc')
parser.add_argument('-r', '--regions', type = str, help= 'BED file containing regions of interest.')
parser.add_argument('-f', '--full', default = True, help="If set to true then count tags for gene/peak")
parser.add_argument('-m', '--force', action='store_true', default=False, help="If set to true then force to count tags again and generate npz")
parser.add_argument('-b', '--bigwig', action='store_true', default=False, help="If set to true then files passed are bigwig files")

args = parser.parse_args()

if args.bigwig is False and args.ctrl is None:
    parser.error("You have not provided a control file and the -b flag is not set.")

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

EXCLUDED_CHROMS = ['chrU', 'chrUextra']

def chrom_filter(f):
    if f.chrom not in EXCLUDED_CHROMS and f.strand != '.':
        return True

def genes_generator(database):# {{{
    """
    Generator function to yield TSS +/- Kb of each annotated gene
    """
    for gene in database.features_of_type('gene'):
        if (re.match('chr', gene.chrom)):
            yield asinterval(gene)

# }}}

# Function Not used -- check if ever needed before trashing
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
            yield TSS(asinterval(transcript), upstream = args.upstream, downstream = args.downstream)
            
# }}}

# For comparing versions.
def cmp_version(version1, version2):
    def normalize(v):
        return [int(x) for x in re.sub(r'(\.0+)*$','', v).split(".")]
    return cmp(normalize(version1), normalize(version2))

# when using filter from pybedtools it does not correctly write the gtf and it can not find the attribute gene_id afterwards# {{{
def correct_bedtools_output(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    f.close()
    f = open(filename, 'w')

    for line in lines:
        f.write(line.strip("\n")+";\n")
    f.close()
# }}}

# Create arrays in parallel, and save to disk for later # {{{
def count_tags(features, description, bw = None, ip = None, ctrl = None):
    "This counts mapped reads for ip and input and normalizes them by library size and million mapped reads"
    assert bw is not None or (ip is not None and ctrl is not None), \
            "Provide either one normalised bigwig or BAM files for IP and control."

    from metaseq import persistence
    import multiprocessing
    # multiprocessing.cpu_count() gives the allocated number of cores, so if used with LSF this will 
    # return the number of cores on the host -- not good practice
    # processes = multiprocessing.cpu_count()
    processes = int(os.environ["LSB_DJOB_NUMPROC"])
    #processes = int(2)
    print 'Counting tags...'

    if bw is not None:
        # Use prefix of ip file as filename for numpy array (npz)
        basename = args.out_dir + os.path.splitext(os.path.basename(bw))[0]
    else:
        basename = args.out_dir + os.path.splitext(os.path.basename(ip))[0]
    # description is a description of the features provided
    basename = basename + '.' + description
    output = basename + '.npz'
    print 'Output: %s' % output
    
    # Calculate number of bins depending on feature# {{{
    pattern = re.compile('(\d+)-(\w+)-(\d+)')
    # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
    matching = pattern.match(description)

    if matching:
        upstream, feature_type, downstream = matching.groups()
        # Set the bins such that counting is done for every 10bp window
        bins = (int(upstream) + int(downstream))/10
    elif re.compile('(\d+)-(\w+)').match(description) or re.compile('(\w+)-(\d+)').match(description):
        if re.compile('(\d+)-(\w+)').match(description):
            upstream, feature_type = re.compile('(\d+)-(\w+)').match(description).groups()
            bins = int(upstream)/10
        else:
            feature_type, downstream = re.compile('(\w+)-(\d+)').match(description).groups()
            bins = int(downstream)/10
    else:
        # For genes it only makes sense to count every 1% of gene
#        bins = 100
        bins = 1000

   # the division makes the number a float which doesn't work with 
    bins = int(bins)
    # }}}

    # Run if file does not exist and experiment has no replicates
    if not os.path.exists(output) or args.force:
        if bw is not None:
            bw = metaseq.genomic_signal(bw, 'bigwig')
            print "Building array based on provided bigwig file for %s using %s processors" % (basename, processes)
            bw_array = bw.array(features, bins = bins, processes = processes)
            persistence.save_features_and_arrays(
                    features = features,
                    arrays={'bw': bw_array},
                    prefix = basename,
                    link_features = True,
                    overwrite = True)
            print "Done"
        else:
            # Read in bam files
            ip = metaseq.genomic_signal(ip, 'bam')
            ctrl = metaseq.genomic_signal(ctrl, 'bam')
            
            # Create arrays in parallel
            print "Building the IP array for %s using %s processors" % (basename, processes)
            ip_array = ip.array(features, bins = bins, processes = processes, shift_width = int(args.fs)/2)
            print "Building the input array for %s using %s processors" % (basename, processes)
            ctrl_array = ctrl.array(features, bins = bins, processes = processes, shift_width = int(args.fs)/2)
            
            # Normalize to library size
            # operation on the right side is performed first
            ip_array /= ip.mapped_read_count() / 1e6
            ctrl_array /= ctrl.mapped_read_count() / 1e6
            # RPKM would be ip_array /= ip.mapped_read_count() / 1e6 => ip_array *= binsize
            ip_array /= 10
            ctrl_array /= 10

            # Cache to disk (will be saved as "example.npz" and "example.features")
            persistence.save_features_and_arrays(
                    features = features,
                    arrays={'ip': ip_array, 'input': ctrl_array},
                    prefix = basename,
                    link_features = True,
                    overwrite = True)
            print "Done"
    else:
        print ".npz already exists for %s using that" % output
    return;
# }}}

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
        print "Creating genes: %s " % filename
        genes = pybedtools.BedTool(genes_generator(db)).saveas(filename)
    else:
        genes = pybedtools.BedTool(filename)
                
    suffix = '.metaseq.' + str(args.upstream) + '-gene_start-' + str(args.downstream) + '.gtf'
    filename = gff_filename.replace('.gtf', suffix)
    if not os.path.exists(filename):
        # A BedTool made out of a generator, and saved to file.
        print "Creating gene_start: %s " % filename
        gene_start = pybedtools.BedTool(gene_start_generator(db))\
                .saveas(filename)
    else:
        gene_start = pybedtools.BedTool(filename)
    
    filename = gff_filename.replace('.gtf', '.metaseq.2000-gene_start.gtf')
    if not os.path.exists(filename) or os.stat(filename).st_size == 0:
        print 'Creating ' + filename + ' ...'
        print "Creating upstream region of gene start: %s " % filename
        upstream = genes.flank(genome=args.assembly, s=True, r=0, l=2000).saveas(filename)
    else:
        upstream = pybedtools.BedTool(filename)
    
    filename = gff_filename.replace('.gtf', '.metaseq.gene_end-500.gtf')
    if not os.path.exists(filename) or os.stat(filename).st_size == 0:
        print "Creating downstream region of gene end: %s " % filename
        downstream = genes.flank(genome=args.assembly, s=True, l=0, r=500).saveas(filename)
    else:
        downstream = pybedtools.BedTool(filename)
   
    # upstream.filter(greater_than, 4999) return a subset of the bedtool object (only features of length bigger than 4999 are included)
    upstream_filtered = [ g.name for g in upstream.filter(greater_than, 1999) ]
    downstream_filtered = [ g.name for g in downstream.filter(greater_than, 499) ]
    print "Genes: %s Pass filter, upstream: %s downstream: %s" % (len(genes), len(upstream_filtered), len(downstream_filtered))

    tmp = [g.name for g in genes if g.name not in upstream_filtered]
    # Joining the 2 tuples. tmp contains the genes to be filtered out
    tmp = tmp + [g.name for g in genes if g.name not in downstream_filtered]

    if len(tmp) != 0 and not os.path.exists(gff_filename.replace('.gtf', '.metaseq.genes-filtered-out.gtf')):
        print "Filtering out regions of width < 5Kb."
        print "These are for genes near chromosome ends that promoter or downstream region hang over chromosome end. (Total: %s Filtered out: %s)." % (len(genes), len(tmp))
        print "Writing filtered files after matching genes and their order."

        # Keeping track of what was excluded
        filtered_out_genes = genes.filter(lambda gene: gene.name in tmp).saveas(gff_filename.replace('.gtf', '.metaseq.genes-filtered-out.gtf'))

        filtered_genes = genes.filter(lambda gene: gene.name not in tmp).saveas(gff_filename.replace('.gtf', '.metaseq.genes.filtered.gtf'))
        if cmp_version(pybedtools.__version__, '0.6.9') <= 0:
            correct_bedtools_output(gff_filename.replace('.gtf', '.metaseq.genes.filtered.gtf'))

        filtered_upstream = upstream.filter(lambda gene: gene.name not in tmp).saveas(gff_filename.replace('.gtf', '.metaseq.2000-gene_start.filtered.gtf'))
        if cmp_version(pybedtools.__version__, '0.6.9') <= 0:
            correct_bedtools_output(gff_filename.replace('.gtf', '.metaseq.2000-gene_start.filtered.gtf'))

        filtered_downstream = downstream.filter(lambda gene: gene.name not in tmp).saveas(gff_filename.replace('.gtf', '.metaseq.gene_end-500.filtered.gtf'))
        if cmp_version(pybedtools.__version__, '0.6.9') <= 0:
            correct_bedtools_output(gff_filename.replace('.gtf', '.metaseq.gene_end-500.filtered.gtf'))
    else:
        filtered_genes = pybedtools.BedTool(gff_filename.replace('.gtf', '.metaseq.genes.filtered.gtf'))
        filtered_upstream = pybedtools.BedTool(gff_filename.replace('.gtf', '.metaseq.2000-gene_start.filtered.gtf'))
        filtered_downstream = pybedtools.BedTool(gff_filename.replace('.gtf', '.metaseq.gene_end-500.filtered.gtf'))
        
    return {'tss':tss, 'gene_start':gene_start, 'genes':genes, 'genes-filtered':filtered_genes, 'upstream-filtered':filtered_upstream, 'downstream-filtered':filtered_downstream}
# }}}

# }}}

# For each bam calculate signal # {{{
def main():

    print 'Calling create_features() ...'
    features = create_features()
    descriptions = {'gene_start': str(args.upstream) + '-gene_start-' + str(args.downstream),
                    'tss': str(args.upstream) + '-tss-' + str(args.downstream),
                    'genes':'genes',
                    'genes-filtered':'genes-filtered', 
                    'upstream-filtered': '2000-gene_start', 
                    'downstream-filtered':'gene_end-500'}

    # Create genomic_signal objects that point to data files
    for i in range(len(args.ip)):

        if args.bigwig:
            for key, value in descriptions.iteritems():
                count_tags(bw = args.ip[i], features = features[key], description = value)
        else:
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

            for key, value in descriptions.iteritems():
                count_tags(ip = bams['ip'], ctrl = bams['ctrl'],\
                        features = features[key], description = value)

        if not args.regions:
            print 'Regions not provided'
        else:
            'ToDo: read in the peaks and do similar analysis to genes'

# }}}

main()
