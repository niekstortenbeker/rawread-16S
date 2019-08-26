from pathlib import Path


basepath = Path(__file__).resolve().parent.parent
output_path = basepath / 'results'

def main():
    files = get_fasta_files(output_path)
    remove_space_fasta(files, '_n')

def get_fasta_files(path):
    return list(path.glob('*.fasta'))

def remove_space_fasta(paths, suffix):
    for path in paths:
        out_path = Path(path.parent, path.stem + suffix + path.suffix)
        with open(path, 'r') as fin, open(out_path, 'w') as fout:
            for line in fin:
                if line[0] == '>':
                    fout.write(line.replace(' ', ''))
                else:
                    fout.write(line)

if __name__ == '__main__':
    main()
