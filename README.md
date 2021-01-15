# EPSim


EP-sim: Multiple-resolution alignment-free measure based on Entropic Profiles

Abstract:

The use of fast similarity measure, like alignment-free measures, to detect related regulatory sequences is crucial to understand functional correlation between two enhancers. However alignment-free measures are generally tied to a fixed resolution $k$. Here we propose an alignment-free statistic, called EP2*, that is based on multiple resolution patterns derived from Entropic Profiles. Entropic Profile is a function of the genomic location that captures the importance of that region with respect to the whole genome. We evaluate several alignment-free statistics on simulated data and real mouse ChIP-seq sequences. The new statistic, EP2*, is highly successful in discriminating functionally related enhancers and, in almost all experiments, it outperforms fixed-resolution methods.

Software

Here you can find the C++ application EP-sim with some examples.

EP-sim is based on the library SEQAN

Run EP-sim using the command:

ep_sim FastaFile OutputFile -D D2* -L k -P 7 -B 1

Using the input sequences contained in "FastaFile" EP-sim computes the statistics EP*2 for all pairs of sequences in the input file, and outputs the results in "OutputFile" as matrix of scores. The length of patterns used is up to k, with variance 0.7 and a background Markov model of order 1.

To run the example included type:

./ep_sim EbolaMixedGuinea.fasta ebolamix_l5_b1_EP2s.txt -D D2* -L 5 -P 7 -B 1 -F -R > junk.txt


Licence

The software is freely available for academic use.
For questions about the tool, please contact Matteo Comin or Morris Antonello .


Reference

Please cite the following papers:

M. Comin, M. Antonello
"On the comparison of regulatory sequences with multiple resolution Entropic Profiles"
Accepted in BMC Bioinformatics, 2016. Pdf

 

 
