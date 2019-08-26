from pathlib import Path
import subprocess as sp
import sys


"""should be run as 'python sortmeRNA.py nslots wdir'. The script should be
in a script directory, the parent directory of which should contain a directory
'data' (can by a symbolic link) output from readQC_rawread.py."""


nslots = sys.argv[1]
wdir = Path(sys.argv[2])

basepath = Path(__file__).resolve().parent.parent
data_path = basepath / 'data'
output_path = basepath / 'results'
output_path.mkdir(exist_ok=True)

print('make temp directories')
tempdir = wdir / 'temp'
tempdir.mkdir()
resultsdir = wdir / 'results'
resultsdir.mkdir()

def main():
    """ functions return paths as Path objects"""
    print('run main')
    data_directories = get_data_directories(data_path)
    for sequence_file in data_directories:
        copied = copy(sequence_file)
        gunzipped = gunzip(copied)
        sortmerna(gunzipped)
        # TODO output/convert fastq to fasta file
        # TODO remove spaces in fasta file 
        # TODO gzip *blast and *fastq
        copy_all_files(resultsdir, output_path)


def get_data_directories(sequence_path):
    """get a list with the directories to the sequence files in sequence_path."""
    return list(data_path.glob('*.fastq.gz'))


def copy(file):
    outdir = tempdir / file.name
    sp.run(['cp', get_path(file), get_path(outdir)])
    return outdir


def gunzip(file):
    """note gzip removes the original file"""
    sp.run(['gzip',
            '-d',
            get_path(file)])
    return file.with_suffix('')  # remove the last .gz suffix


def sortmerna(reads_file):
    """note this will only work if the shell script runs "module load sortmerna
    Hmm. Oh no that is probably not true"""
    dbs = '/bioinf/software/sortmerna/sortmerna-2.0/rRNA_databases'
    index = '/bioinf/software/sortmerna/sortmerna-2.0/index'
    databases = f'{dbs}/silva-bac-16s-id90.fasta,{index}/silva-bac-16s-db:{dbs}/silva-arc-16s-id95.fasta,{index}/silva-arc-16s-db:{dbs}/silva-euk-18s-id95.fasta,{index}/silva-euk-18s-db'
    out_base = resultsdir / f'{reads_file.stem}_rRNA'
    logdir = resultsdir / f'{reads_file.stem}_sortmeRNA_log.txt'
    # TODO I don't really want those blast results, can I also not output them
    with open(logdir, 'w') as logfile:
        sp.run(['sortmerna',
                '--ref', databases,
                '--reads', get_path(reads_file),
                '--aligned', get_path(out_base),
                '--fastx',  # output fasta/fastq file
                '--blast', '1',  # output alignments in BLAST tabular format
                '--log',
                '--num_alignments', '0',  #all alignments will be output
                '-v',  # verbose
                '-a', nslots],
                stdout=logfile, stderr=sp.STDOUT)


def copy_all_files(fin, fout):
    """accepts either Pathlib paths or str. If is not str, assumed to be a
    pathlib.Posixpath"""
    fin = get_path(fin)
    fout = get_path(fout)
    sp.run(f'cp {fin}/* {fout}', shell=True)


def get_path(path):
    """get path as a str. If not already a str it is assumed to be
    a pathlib.Path object"""
    if isinstance(path, str):
        return path
    else:
        return str(path.resolve())


if __name__ == '__main__':
    main()
