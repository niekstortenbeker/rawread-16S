# rawread-16S
Shell and python scripts to automatically analyze 16S sequences from raw illumina sequencing files.
The scripts automatically recognize libraries etc. based on the file names that are generated by 
illumina sequencing. 

The python scripts do all the file handeling and execution of command line programs. The shell scripts 
execute the python scripts on a cluster. 
First the readQC_rawread should run, and then the sortmerna. When necessary fasta headers can be cleaned
with remove_space_fasta.py. This data can then be uploaded to silvaNGS. Scripts results and data should 
be in a file structure with separate directories for scripts, data and results like so:

-readQC_rawread
----------------scripts
----------------results
----------------data
-sortmerna
----------------scripts
----------------results
----------------data
