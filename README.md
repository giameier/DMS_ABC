# Deep mutational scanning of ABC transporter

Deep mutational scanning of ABCtransporter is a software for processing and analyzing of sequence alignment files from deep mutational scanning experiments of ABC transporter unsing paired overlapp sequences.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites and Installing

Enrich2 runs on Python 3.7 and needs the following dependencies:
* pysam version 0.15.3
* Tkinter


To install create a environment and activate it

download the code

Execute the main.py script in your environment which starts the graphical user interface


## Running 
To start a new project the following parameters have to be set:
* Input file directory: This directory should contain all your input files to be analyzed. (supported are bam files with endings ".sorted.bam" and the folder should also contain an index file ".sorted.bam.bai")
* Output file dierctory: Directory where the output files should be saved
* Job Name: A name for your project
* Positions on the gene to be analyzed: This should be given as the position of the first base of codon. Positions should be comaseparated and contain no spaces
* DNA sequence of the reference: This is your gene of interrest as it was used during alingment with BWA. (Only ATGC allowed and no lower case)
* Position of last vase of reading frame1: If you have several genes on your reference sequence you have to specify the last base of the first frame
* Frameshift offset: frame shift of secondond frame relative to first frame.
* If you like to load parameters from an old project, specify the location of the respective Json file 
The programm supports parallelprocessing of files. You can specify the amount of CPUs to be used in the DMS_processing_multiprocessing.py file.


After hitting the process button the following is executed:

All sorte.bam files in the input directory are found.
All parameters are saved to a Json file in the output directory.
For each input file a bam file is created containing trimmed sequences
For each input file a variant count file and HDF5 files are created

## Authors

**Gianmarco Meier** University of Zurich 

