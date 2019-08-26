from pathlib import Path
import re
import subprocess as sp
import sys

"""should be run as 'python sortmeRNA.py nslots wdir'. The script should be
in a script directory, the parent directory of which should contain a directory
'data' (can by a symbolic link) output from readQX_rawread.py. This script will
output a fastq.gz file with the reads ready for sortmeRNA (or comparable stuff),
and FASTQC reports"""


nslots = sys.argv[1]
wdir = Path(sys.argv[2])

basepath = Path(__file__).resolve().parent.parent
sequence_path = basepath / 'data'
sequence_path = sequence_path.resolve()
print(sequence_path)
output_path = basepath / 'results'
output_path.mkdir(exist_ok=True)
print(output_path)

print('make temp directories')
tempdir = wdir / 'temp'
tempdir.mkdir()
resultsdir = wdir / 'results'
resultsdir.mkdir()


def main():
    """every time the functions return paths of the files that are created in
    the function. Note, these paths are as strings, not path objects"""
    print('Okay, start running the script')
    sequence_dictionary, libraries = get_sequence_directories(sequence_path)
    for library in libraries:
        zcatted_files = zcatting_runs(library, sequence_dictionary)
        trimmed_reads = trim_reads(library, zcatted_files, nslots)
        merged_reads = merge_reads(library, trimmed_reads, nslots)
        appended_reads = append_reads(library, merged_reads)
        fastqc(appended_reads)
        gzip(appended_reads)
        copy_all_files(str(resultsdir.resolve()), output_path)
        empty_directories([resultsdir, tempdir])


def get_sequence_directories(sequence_path):
    """get a dictionary that saves the directories to the sequence files.
    Keys are library names with a fw or rv suffix like so: '4111_A_R1'. The
    values are lists with one or more runs beloning to each library and sequence
    direction. The libraries is a list with all the library names"""
    sequence_files = list(sequence_path.glob('*.fastq.gz'))
    libraries = set([sequence_file.name[0:6] for sequence_file in sequence_files]) # e.g. 4111_A
    fw_rv = ['_R1_001', '_R2_001']

    sequence_dictionary = {}
    for library in libraries:
        for direction_suffix in fw_rv:
            runs = extract_illumina_run(library, direction_suffix, sequence_files)
            ID = library + '_' + direction_suffix[1:3]  #e.g. 4111_ + _ + R1
            sequence_dictionary[ID] = runs
    return sequence_dictionary, libraries


def extract_illumina_run(run_prefix, direction_suffix, sequence_files):
    runs = []
    for sequence_file in sequence_files:
        file_name = sequence_file.name
        match_str = f'{run_prefix}(.*){direction_suffix}(.*)'
        match = re.match(match_str, file_name)
        if match:
            runs.append(str(sequence_file.resolve()))
    return runs


def zcatting_runs(library, sequence_dictionary):
    """zcat sequence runs. When a library consists of multiple runs (as is
    determined in the sequence_dictionary) concanate these files. returns
    a list with [fw_path, rv_path] (e.g. path/4441_A_R1.fastq) """
    output = []
    for direction in ['_R1', '_R2']:
        files = ' '.join(sequence_dictionary[library+direction])
        new_file = library + direction + '.fastq'
        new_file = tempdir / new_file
        new_file = str(new_file.resolve())
        output.append(new_file)
        sp.run(f'zcat {files} > {new_file}', shell=True)
    return output


def trim_reads(library, zcatted_files, nslots):
    """QC: trim reads with low quality base calling. Return list with paths as
    strings for merged, fw, rv."""
    fw = zcatted_files[0]
    rv = zcatted_files[1]
    base = f'{str(tempdir.resolve())}/{library}'
    out_fw = f'{base}_fw_bbduk.fastq'
    out_rv = f'{base}_rv_bbduk.fastq'
    base = f'{str(resultsdir.resolve())}/{library}'
    logdir_paired = f'{base}_bbduk_log.txt'
    BBDuk_paired(fw, rv, out_fw, out_rv, logdir_paired, nslots)
    output = [out_fw, out_rv]
    return output


def BBDuk_paired(fwin, rvin, fwout, rvout, logdir, nslots):
    """QC of the paired end reads. trim read ends that are low quality and general
    read areas below a quality threshold. returns the trimmed reads directory as
    strings)"""
    BBDuk_path = '/bioinf/software/bbmap/bbmap-35.14/bbduk.sh'
    with open(logdir, 'w') as logfile:
        sp.run([BBDuk_path,
                f'in={fwin}',
                f'in2={rvin}',
                f'out={fwout}',
                f'out2={rvout}',
                'qin=33',  # input quality offset
                'minlen=30',  # Reads shorter than this after trimming will be discarded.
                'qtrim=rl',  # remove low quality ends, trim both ends
                'trimq=20',  # Regions with average quality BELOW this will be trimmed. (q20 is 1/100 is false)
                f't={nslots}'],  # no threads
                stdout=logfile, stderr=sp.STDOUT)


def merge_reads(library, trimmed_reads, nslots):
    """merge reads with PEAR. Return the output files,
    as a list with merged, fw and rv reads directories as strings"""
    fw = trimmed_reads[0]
    rv = trimmed_reads[1]
    output_base = f'{str(tempdir.resolve())}/{library}'
    logdir = f'{str(resultsdir.resolve())}/{library}_PEAR_log.txt'
    PEAR(fw, rv, output_base, logdir, nslots)
    output = [f'{output_base}.assembled.fastq',
              f'{output_base}.unassembled.forward.fastq',
              f'{output_base}.unassembled.reverse.fastq']
    return output


def PEAR(forward, reverse, output_base, logdir, nslots):
    """merge reads. uses standard settings. p=0.01, overlap -v = 10 bp
    outputs {output_base}assembled.fastq {output_base}unassembled.forward.fastq
    {output_base}unassembled.reverse.fastq {output_base}discarded.fastq"""
    PEAR_path = '/bioinf/software/pear/pear-0.9.8/bin/pear'
    with open(logdir, 'w') as fout:
        sp.run([PEAR_path,
                '-f', forward,
                '-r', reverse,
                '-o', output_base,  # Specify the name to be used as base for the output files (can be with directory). PEAR outputs four files. A file containing the assembled reads with a assembled.fastq extension, two files containing the forward, resp. reverse, unassembled reads with extensions unassembled.forward.fastq, resp. unassembled.reverse.fastq, and a file containing the discarded reads with a discarded.fastq extension.
                '-j', nslots],  # no threads
                stdout=fout, stderr=sp.STDOUT)


def append_reads(library, reads):
    """use cat to append reads. Returns a string with the
    output directory"""
    fin = ' '.join(reads)
    fout = str(resultsdir.resolve())
    fout = f'{fout}/{library}.fastq'
    sp.run(f'cat {fin} > {fout}', shell=True)
    return fout


def fastqc(fastq_file):
    """generate fastqc file for one fastq file, in the same directory as the
    input file"""
    fastqc = '/bioinf/software/fastqc/fastqc-0.11.4/fastqc'
    sp.run([fastqc,
            '-t', nslots,
            '--noextract',  # do not decompress output files
            '--quiet',
            fastq_file])


def gzip(file):
    """note gzip removes the original file"""
    sp.run(['gzip',
             file])


def get_path(path):
    """get path as a str. If not already a str it is assumed to be
    a pathlib.Path object"""
    if isinstance(path, str):
        return path
    else:
        return str(path.resolve())


def copy_all_files(fin, fout):
    """accepts either Pathlib paths or str. If is not str, assumed to be a
    pathlib.Posixpath"""
    fin = get_path(fin)
    fout = get_path(fout)
    sp.run(f'cp {fin}/* {fout}', shell=True)


def empty_directories(directories):
    directories = [f'{get_path(directory)}/*' for directory in directories]
    directories = ' '.join(directories)
    sp.run(f'rm {directories}', shell=True)


if __name__ == '__main__':
    main()
