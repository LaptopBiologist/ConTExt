# ConTExt
A bioinformatic pipeline for identifying structural variation in repetitive DNA using NGS data.

Repeat-derived NGS reads are often ignored in sequencing analyses because an essential step in most NGS pipelines is identifying the genomic locus from which a given read originates. For reads originating in repetitive sequence, this is often impossible. Consequently, asking questions about repeats using NGS data requires phrasing the questions differently than one would for unique sequence. ConTExt is a set of tools designed to reframe questions in this manner. Instead of trying to identify the locus from which a repeat-derived read originates, it seeks to determine the repeat-family from which the read originates and the location within the repeat-family (e.g. for a TE-derived read, does it sit near the 5' end of the element or the nearer to the 3' end?). Once reads are organized in this manner, ConTExt makes use of read pair information, coverage, and sequence polymorphism to identify structures involving repetitive sequence (tandem junctions, deletions within a repeat, insertions into unique or repeated sequence), estimate their copy number, and identify the proportion repeats harboring variant alleles. 

Most of the major scripts that actually organize the data are designed to be run from command-prompt and all of the necessary steps are called automatically. The individual functions are not intended to be called directly by the user. Extensive documentation for using functions is restricted to those intended to be called manually.  

I automate as much of the pipeline as is possible, but cookie cutter solutions do always not exist for every step of analysis. Where user decision are necessary, the manual provided in ./docs provides guidance. 

Dependencies:
Python 2.7, the code is not currently compatible with Python 3. In future, I may try to port it to Python 3, but it is not currently a priority of mine. The package was developed using the following versions of these Python libraries:

numpy>=1.10.4, scipy>=0.17.1, Biopython>=1.6.8, matplotlib>=1.5.1, scikit-learn>=0.18.1, seaborn>=0.7.1

It may run with earlier versions of some of these, but not all. If you have earlier versions installed, the provided setup.py script will update your packages.


Installation:
Download the repository and run the setup.py script, which will check for and install the requisite dependencies:

(in the directory containing setup.py)

python setup.py install


