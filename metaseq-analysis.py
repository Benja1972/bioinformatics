#!/usr/bin/env python

#~ https://pythonhosted.org/metaseq/example_session.html

import os
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval


ref_dir = '../references'
data_dir = '04aln'



db = gffutils.FeatureDB(os.path.join(ref_dir,'mm10_UCSC_EnsemblGenes_ensGene.gtf.db'))


def tss_generator():
    """
    Generator function to yield TSS of each annotated transcript
    """
    for transcript in db.features_of_type('transcript'):
        yield TSS(asinterval(transcript), upstream=1, downstream=0)


# A BedTool made out of a generator, and saved to file.
tsses = pybedtools.BedTool(tss_generator()).saveas(os.path.join(ref_dir,'tsses.gtf'))

tsses_1kb = tsses.slop(b=1000, genome='mm10', output=os.path.join(ref_dir,'tsses-1kb.gtf'))

# ---------------------------------------------------------
#~ Metaseq #
# ---------------------------------------------------------

#~ metaseq works with the concepts of signal and windows. In this example, the signal is ChIP data, and the windows are TSS +/- 1kb.

#~ The first step is to create “genomic signal” objects out of the data. Since our example files are BAM files, we specify the kind=’bam’, but if you have your own data in a different format (bigWig, bigBed, BED, GFF, GTF, VCF) then specify that format instead (see metaseq.genomic_signal()).

import metaseq

ip_file = 'TLX3_H3K4me3_repl1.sorted.bam'
#ip_file = 'TLX3_H3K4me3_repl2.sorted.bam'


input_file = 'INP-TLX3_1.sorted.bam'
#input_file = 'INP-TLX3_2.sorted.bam'


ip_signal = metaseq.genomic_signal(os.path.join(data_dir, ip_file),'bam')

input_signal = metaseq.genomic_signal(os.path.join(data_dir, input_file),'bam')

#~ Now we can create the arrays of signal over each window. Since this can be a time-consuming step, the first time this code is run it will cache the arrays on disk. The next time this code is run, it will be quickly loaded. Trigger a re-run by deleting the .npz file.

#~ Here, with the BamSignal.array method, we bin each promoter region into 100 bins, and calculate the signal in parallel across as many CPUs as are available. We do this for the IP signal and input signals separately. Then, since these are BAM files of mapped reads, we scale the arrays to the library size. The scaled arrays are then saved to disk, along with the windows that were used to create them.


import multiprocessing
processes = multiprocessing.cpu_count()

if not os.path.exists('example.npz'):

    # The signal is the IP ChIP-seq BAM file.
    ip_array = ip_signal.array(

        # Look at signal over these windows
        tsses_1kb,

        # Bin signal into this many bins per window
        bins=100,

        # Use multiple CPUs. Dramatically speeds up run time.
        processes=processes)

    # Do the same thing for input.
    input_array = input_signal.array(
        tsses_1kb,
        bins=100,
        processes=processes)

    # Normalize to library size. The values in the array
    # will be in units of "reads per million mapped reads"
    ip_array /= ip_signal.mapped_read_count() / 1e6
    input_array /= input_signal.mapped_read_count() / 1e6

    # Cache to disk. The data will be saved as "example.npz" and "example.features".
    metaseq.persistence.save_features_and_arrays(
        features=tsses,
        arrays={'ip': ip_array, 'input': input_array},
        prefix='example',
        link_features=False,
        overwrite=True)

features, arrays = metaseq.persistence.load_features_and_arrays(prefix='example')


# How many features?
#assert len(features) == 5708

# This ought to be exactly the same as the number of features in `tsses_1kb.gtf`
assert len(features) == len(tsses_1kb) 

# This shows that `arrays` acts like a dictionary
assert sorted(arrays.keys()) == ['input', 'ip']

# This shows that the IP and input arrays have one row per feature, and one column per bin
assert arrays['ip'].shape == arrays['input'].shape

# ---------------------------------------------------------
#~ Line plot of average signal
# ---------------------------------------------------------

#~ Now that we have NumPy arrays of signal over windows, there’s a lot we can do. One easy thing is to simply plot the mean signal of IP and of input. Let’s construct meaningful values for the x-axis, from -1000 to +1000 over 100 bins. We’ll do this with a NumPy array.

import numpy as np
x = np.linspace(-1000, 1000, 100)

# Import plotting tools
from matplotlib import pyplot as plt


# Create a figure and axes
fig = plt.figure()
ax = fig.add_subplot(111)


# Plot the IP:
ax.plot(
    # use the x-axis values we created
    x,

    # axis=0 takes the column-wise mean, so with
    # 100 columns we'll have 100 means to plot
    arrays['ip'].mean(axis=0),

    # Make it red
    color='r',

    # Label to show up in legend
    label='IP')


# Do the same thing with the input
ax.plot(
    x,
    arrays['input'].mean(axis=0),
    color='k',
    label='input')


# Add a vertical line at the TSS, at position 0
ax.axvline(0, linestyle=':', color='k')


# Add labels and legend
ax.set_xlabel('Distance from TSS (bp)')
ax.set_ylabel('Average read coverage (per million mapped reads)')
ax.legend(loc='best');


# ---------------------------------------------------------
#~ Adding a heatmap
# ---------------------------------------------------------

#~ We don’t really know if this average signal is due to a handful of really strong peaks, or if it’s moderate signal over many peaks. So one improvement would be to include a heatmap of the signal over all the TSSs.

#~ First, let’s create a single normalized array by subtracting input from IP:


normalized_subtracted = arrays['ip'] - arrays['input']

#~ metaseq comes with some helper functions to simplify this kind of plotting. The metaseq.plotutils.imshow function is one of these; here the arguments are described:

# Tweak some font settings so the results look nicer
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# the metaseq.plotutils.imshow function does a lot of work,
# we just have to give it the right arguments:
fig1 = metaseq.plotutils.imshow(

    # The array to plot; here, we've subtracted input from IP.
    normalized_subtracted,

    # X-axis to use
    x=x,

    # Change the default figure size to something smaller for this example
    figsize=(3, 7),

    # Make the colorbar limits go from 5th to 99th percentile.
    # `percentile=True` means treat vmin/vmax as percentiles rather than
    # actual values.
    percentile=True,
    vmin=5,
    vmax=99,

    # Style for the average line plot (black line)
    line_kwargs=dict(color='k', label='All'),

    # Style for the +/- 95% CI band surrounding the
    # average line (transparent black)
    fill_kwargs=dict(color='k', alpha=0.3),
)

# ---------------------------------------------------------
#~ Sorting the array
# ---------------------------------------------------------

#~ The array is not very meaningful as currently sorted. We can adjust the sorting this either by re-ordering the array before plotting, or using the sort_by kwarg when calling metaseq.plotutils.imshow. Let’s sort the rows by their mean value:

fig2 = metaseq.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=normalized_subtracted.mean(axis=1)
)
# "line_axes" is our handle for working on the lower axes.
# Add some nicer labels.
fig2.line_axes.set_ylabel('Average enrichment');
fig2.line_axes.set_xlabel('Distance from TSS (bp)');

# "array_axes" is our handle for working on the upper array axes.
# Add a nicer axis label
fig2.array_axes.set_ylabel('Transcripts on chr17')

# Remove the x tick labels, since they're redundant
# with the line axes
fig2.array_axes.set_xticklabels([])

# Add a vertical line to indicate zero in both the array axes
# and the line axes
fig2.array_axes.axvline(0, linestyle=':', color='k')
fig2.line_axes.axvline(0, linestyle=':', color='k')

fig2.cax.set_ylabel("Enrichment")



#~ We can use any number of arbitrary sorting methods. For example, this sorts the rows by the position of the highest signal in the row. Note that the line plot, which is the column-wise average, remains unchanged since we’re still using the same data. The rows are just sorted differently.

fig3 = metaseq.plotutils.imshow(

    # These are the same arguments as above.
    normalized_subtracted,
    x=x,
    figsize=(3, 7),
    vmin=5, vmax=99,  percentile=True,
    line_kwargs=dict(color='k', label='All'),
    fill_kwargs=dict(color='k', alpha=0.3),

    # This is new: sort by mean signal
    sort_by=np.argmax(normalized_subtracted, axis=1)
)

# ---------------------------------------------------------
#~ RNA-seq data wrangling: loading data
# ---------------------------------------------------------

#~ The metaseq.results_table module has tools for working with this kind of data (for example, the metaseq.results_table.DESeq2Results class). Here, we will make a generic ResultsTable which handles any kind of tab-delimited data. It’s important to specify the index column. This is the column that contains the transcript IDs in these files.


#~ rna_data_dir = '/home/sergio/media/NAS3B/Bioinfoplatform/PFlab/RNA-Seq/RNA-SEQ_ANALYSIS/3-HTSCOUNT/TopHat_Very-Sensitive'

#~ control_file = 'RAGZ.HTSeqCount.txt'
#~ knockdown_file = 'TLX3-1_1.HTSeqCount.txt'


rna_data_dir = '/home/sergio/media/NAS3B/Bioinfoplatform/PFlab/RNA-Seq/RNA-SEQ_ANALYSIS/4-DEG_DESeq-edgeR/TopHat_Very-Sensitive/Results'


control_file = 'RAG_VS_TLX3_DEG_By_edgeR.txt'
knockdown_file = 'RAG_VS_TAP_DEG_By_edgeR.txt'





from metaseq.results_table import ResultsTable, DESeq2Results

#~ control = ResultsTable(
    #~ os.path.join(rna_data_dir, control_file),
    #~ import_kwargs=dict(index_col=0, names=['score']))

#~ knockdown = ResultsTable(
    #~ os.path.join(rna_data_dir, knockdown_file),
    #~ import_kwargs=dict(index_col=0, names=['score']))

control = DESeq2Results(
    os.path.join(rna_data_dir, control_file),
    import_kwargs=dict(index_col=0))

knockdown = DESeq2Results(
    os.path.join(rna_data_dir, knockdown_file),
    import_kwargs=dict(index_col=0))


# Inspect results to see what we're working with

print len(control.data)
control.data.head()

# ---------------------------------------------------------
#~ RNA-seq data wrangling: aligning RNA-seq data with ChIP-seq data
# ---------------------------------------------------------

#~ We should ensure that control and knockdown have their transcript IDs in the same order as the rows in the heatmap array, and that they only contain transcript IDs from chr17.

#~ The ResultsTable.reindex_to method is very useful for this – it takes a pybedtools.BedTool object and re-indexes the underlying dataframe so that the order of the dataframe matches the order of the features in the file. In this way we can re-align RNA-seq data to ChIP-seq data for more direct comparison.

#~ Remember the tsses_1kb object that we used to create the array? That defined the order of the rows in the array. We can use that to re-index the dataframes. Let’s look at the first line from that file to see how the transcript ID information is stored:

# Inspect the GTF file originally used to create the array

print tsses_1kb[0]

#~ The ResultsTable is indexed by transcript ID. Note that DESeq2 and edgeR results are typically indexed by gene, rather than trancscript, ID. So when working with your own data, be sure to select the GTF attribute whose values will be found in the ResultsTable index.

#~ Here, we tell the ResultsTable.reindex_to method which attribute it should pay attention to when realigning the data:

# Re-align the ResultsTables to match the GTF file
control = control.reindex_to(tsses, attribute='transcript_id')
knockdown = knockdown.reindex_to(tsses, attribute='transcript_id')

# Everything should be the same length
assert len(control.data) == len(knockdown.data) == len(tsses_1kb)

# Spot-check some values to make sure the GTF file and the DataFrame match up.
assert tsses[0]['transcript_id'] == control.data.index[0]
assert tsses[100]['transcript_id'] == control.data.index[100]
assert tsses[5000]['transcript_id'] == control.data.index[5000]

# ---------------------------------------------------------
#~ RNA-seq data wrangling: join control and knockdown data
# ---------------------------------------------------------

#~ Now for some more data-wrangling. We’ll use basic ``pandas` <http://pandas.pydata.org/>`__ operations to merge the control and knockdown data together into a single table. We’ll also create a new log2foldchange column.

# Join the dataframes and create a new pandas.DataFrame.
data = control.data.join(knockdown.data, lsuffix='_control', rsuffix='_knockdown')

# Add a log2 fold change variable
#~ data['log2foldchange'] = np.log2(data.fpkm_knockdown / data.fpkm_control)
#~ data.head()

# ---------------------------------------------------------
# How many transcripts on chr17 changed expression?

#~ print "up:", sum(data.log2foldchange > 1)
#~ print "down:", sum(data.log2foldchange < -1)
