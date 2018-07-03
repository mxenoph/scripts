
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
parser.add_argument('-r', '--rep', action='store_true', default=False, help="If set to true then pattern marches replicates")
parser.add_argument('-l', '--series', action='store_true', default=False, help="If set to true then pattern marches a series of experiments")
parser.add_argument('-c', '--colors', nargs = '+', type = str, required = False,  help = 'Pattern:color string')

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



# }}}

def plot_metagene():
    array_order = ['2000-gene_start.npz', 'genes-filtered.npz', 'gene_end-500.npz']
    tmp = [ os.path.join(args.path, f) for f in os.listdir(args.path) if f.endswith(".npz") and fnmatch.fnmatch(f, args.matching) ]
    #tmp = [ f for f in tmp if f.endswith('.5000-gene_start.npz') or f.endswith('.genes-filtered.npz') or f.endswith('.gene_end-5000.npz') ]
    tmp = [ f for f in tmp if f.endswith('.2000-gene_start.npz') or f.endswith('.genes-filtered.npz') or f.endswith('.gene_end-500.npz') ]
    assert len(tmp) == 3, \
        'Input arrays provided are not 3.'

# Not the best way but only way I can ensure the order of the arrays in the tuple is gene start, gene, gene end
    files = []
    for m in array_order:
        for f in tmp:
            if f.endswith(m):
                files = files + [f,]


# # For each bam calculate signal and plot it# {{{
# # The windows we'll get signal over
from matplotlib.backends.backend_pdf import PdfPages
normalised_arrays = dict()
# Create genomic_signal objects that point to data files
for f in files:
    npz = os.path.splitext(f)[0]
    # Load the windows and arrays
    from metaseq import persistence
    features, arrays = persistence.load_features_and_arrays(prefix = npz)
   
    if 'bw' in arrays:
        normalised_arrays[os.path.basename(npz)] = arrays['bw']
    else:
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
# {{{
if args.series:
    for x in names:
        npz =  [ k for k in normalised_arrays.keys() if re.match(x,k) ]
        for k in npz:
            n = k
            if re.compile('(\d)+-tss-(\d+)').match(os.path.splitext(k)[1].replace(".","")):
                print 'tss'
                per_feature_type['tss'][n] = normalised_arrays.get(n)
                #per_feature_type['tss'].add(n)
            if re.compile('(\d)+-gene_start-(\d+)').match(os.path.splitext(k)[1].replace(".","")):
                print 'gene-strat'
                per_feature_type['gene_start'][n] = normalised_arrays.get(n)
                #per_feature_type['gene_start'].add(n)
            if os.path.splitext(k)[1].replace(".","") == 'genes':
                print 'gene'
                #per_feature_type['genes'].add(n)
                per_feature_type['genes'][n] = normalised_arrays.get(n)# }}}

if args.rep:
# Per replicate print signal over all features # {{{
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
                bins = (int(upstream) + int(downstream))/10
                # bins are of size 100bp in count-tags.py
                x = np.linspace(-int(upstream), int(downstream), bins)
                #x = np.linspace(-int(upstream), int(downstream), 100)
                
                if feature_type == 'tss':
                    per_feature_type['tss'][n] = normalised_arrays.get(n)
                elif feature_type == 'gene_start':
                    per_feature_type['gene_start'][n] = normalised_arrays.get(n)
            elif window == 'genes':
                feature_type = 'genes'
                # gene array goes from o) to 100% in bins of 100bp
                x = np.linspace(0, 100, 1000)

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
    window = [ w.split('.')[1] for w in per_feature_type.get(key) ]
    print per_feature_type.get(key)

    if len(set(window)) == 1:

        def get_feature_type(string):
            pattern = re.compile('(\d+)-(\w+)-(\d+)')
            # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
            matching = pattern.match(string[0])
            upstream = downstream = Non
                x = np.linspace(-int(upstream), int(downstream), ((int(upstream) + int(downstream))/10))
            
            elif re.compile('(\d+)-(\w+)').match(string[0]) or re.compile('(\w+)-(\d+)').match(string[0]):
                
                if re.compile('(\d+)
                    bins = int(upstream)/10
                
                else:
                    feature_type, downstream = re.compile('(\d+)-(\w+)').match(string[0]).groups()
                    bins = int(downstream)/10
            else:
                # gene array goes from o) to 100% in bins of 100bp
                x = np.linspace(0, 100, 1000)

            return {'x': x, 'feature_ty
            features_from_groups = subsets.index
            subsets = subsets.reindex_to(features, attribute='gene_id')
            cls = np.zeros(len(array)).astype('str')

#            inds = []
#            for cls in subset_order:
#                subset_ind = np.nonzero(subset_by == cls)[0
#                        subset_sort_by = sort_by[subset_ind]
#                                subset_argsort_by = np.argsort(subset_sort_by)
#                                        inds.append(subset_ind[subset_argsort_by])
#                                            ind = np.concatenate(inds)

#            plot_average(per_feature_type[key], x, pp = pp, feature_type = key, subsets = args.subsets)
    else:
        print 'Can not plot experiment together. Tags calculated on different windows'
        
pp.close()
# }}}
#for key in per_feature_type.keys():
#    # key %in% tss, gene_start, genes
#    print key
#    window = [ w.split('.')[1] for w in per_feature_type.get(key).keys() ]
#    print per_feature_type.get(key).keys()
#    if len(set(window)) == 1:
#        pattern = re.compile('(\d+)-(\w+)-(\d+)')
#        # matching will be empty if .genes.features or any other file not in upstream-feature-downstream format
#        matching = pattern.match(window[0])
#        
#        if matching:
#            upstream, feature_type, downstream = matching.groups()
#            # bins are of size 100bp in count-tags.py
#            x = np.linspace(-int(upstream), int(downstream), ((int(upstream) + int(downstream))/10))
#        elif re.compile('(\d+)-(\w+)').match(window[0]) or re.compile('(\w+)-(\d+)').match(window[0]):
#            if re.compile('(\d+)-(\w+)').match(window[0]):
#                upstream, feature_type = re.compile('(\d+)-(\w+)').match(window[0]).groups()
#                bins = int(upstream)/10
#            else:
#                feature_type, downstream = re.compile('(\d+)-(\w+)').match(window[0]).groups()
#                bins = int(downstream)/10
#        else:
#            # gene array goes from o) to 100% in bins of 100bp
#            x = np.linspace(0, 100, 100)
#        plot_average(per_feature_type[key], x, pp = pp, feature_type = key)
#    else:
#        print 'Can not plot experiment together. Tags calculated on different windows'
#        
#pp.close()
#
# }}}

