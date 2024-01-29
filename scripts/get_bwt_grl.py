import os
import glob
import gzip
import subprocess
from Bio import SeqIO
import gzip

##################################################
# sdsl-lite installation - grlbwt uses it
##################################################
'''
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
./install.sh
'''
##################################################
# grlbwt installation
##################################################
'''
git clone git@github.com:ddiazdom/grlBWT.git
cd grlBWT
mkdir build
cd build
cmake ..
make
'''
##################################################

FORMAT_MAPPING = {
    '.fasta': 'fasta',
    '.fa': 'fasta',
    '.fastq': 'fastq',
    '.fq': 'fastq',
    '.fasta.gz': 'fasta',
    '.fa.gz': 'fasta',
    '.fastq.gz': 'fastq',
    '.fq.gz': 'fastq',
}


def is_program_installed(program):
    try:
        # Attempt to execute the program with --version or a similar argument that causes it to exit quickly
        subprocess.run([program, '--version'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except FileNotFoundError:
        return False


def process_sequence_file(file_path, output_file_path, file_format):
    concatenated_sequence = ""
    open_func = gzip.open if file_path.endswith('.gz') else open
    file_format = FORMAT_MAPPING.get(file_format.lower())
    with open_func(file_path, "rt") as file:
        for record in SeqIO.parse(file, file_format):
            sequence = str(record.seq).upper()  # Convert sequence to uppercase
            concatenated_sequence += sequence + '$'

    with open(output_file_path, 'w') as file:
        file.write(concatenated_sequence)


def run_grlbwt(file_path, out_file_path, out_txt_file_path, num_threads):
    # change the -T value to be the basename of the out_file_path
    print(f'grlbwt-cli {file_path} -o {out_file_path} -T {os.path.dirname(out_file_path)} -t {num_threads}')
    os.system(f'grlbwt-cli {file_path} -o {out_file_path} -T {os.path.dirname(out_file_path)} -t {num_threads}')
    #print(f"grl2plain {out_file_path}.rl_bwt {out_txt_file_path}")
    #os.system(f'grl2plain {out_file_path}.rl_bwt {out_txt_file_path}')
    # need to remove one of the extensions from the out_file_path.rl_bwt, if it exists
    # so if the out_file_path is "x_x_x.fastq.gz.txt.rl_bwt", then we want to remove ".txt"
    print(f"grl2plain {out_file_path.rsplit('.', 1)[0]}.rl_bwt {out_txt_file_path}")
    os.system(f"grl2plain {out_file_path.rsplit('.', 1)[0]}.rl_bwt {out_txt_file_path}")


def main(input, num_threads, save_intermediate):
    # input file fastq name ex. "./downloads/x_x_x.fastq.gz"
    in_path_for_fastq_file = input
    extension = ''
    for file_ending in FORMAT_MAPPING:
        if str.endswith(input, file_ending):
            extension = file_ending
            break
    if not extension:
        raise ValueError(f'File must be one of the following formats: {FORMAT_MAPPING.keys()}')

    if not is_program_installed('grlbwt-cli'):
        raise EnvironmentError("grlbwt-cli is not installed. Please install it to continue.")
    if not is_program_installed('grl2plain'):
        raise EnvironmentError("grl2plain is not installed. Please install it to continue.")

    # concatenated input file
    in_path_for_concatenated_string = in_path_for_fastq_file + ".txt"
    # grlbwt processed; auto-adds the .rl_bwt extension. I just need to remove the extension from the input file. Note:
    # there could be many `.` in the file name, so just take the last one.
    #out_path_grlbwt_rl_bwt_file = in_path_for_fastq_file.rsplit('.', 1)[0]
    out_path_grlbwt_rl_bwt_file = in_path_for_fastq_file
    # grlbwt processed to text. ex. ".{file}.bwt"
    out_txt_path_grlbwt_txt_file = in_path_for_fastq_file + ".bwt"
    print(f'input: {input}')
    print(f'sequences will be concatenated : {in_path_for_concatenated_string}')
    print(f'grlbwt will be processed : {out_path_grlbwt_rl_bwt_file}')
    print(f'grlbwt output will be converted to txt : {out_txt_path_grlbwt_txt_file}')

    # preprocess
    if os.path.isfile(in_path_for_concatenated_string):
        print('already concatenated')
    else:
        process_sequence_file(in_path_for_fastq_file, in_path_for_concatenated_string, extension)

    # grlbwt and conversion to txt
    run_grlbwt(in_path_for_concatenated_string, out_path_grlbwt_rl_bwt_file, out_txt_path_grlbwt_txt_file, num_threads)

    # Remove intermediate files
    if not save_intermediate:
        intermediate_files = [in_path_for_concatenated_string, out_path_grlbwt_rl_bwt_file + ".rl_bwt"]
        for file in intermediate_files:
            if os.path.exists(file):
                os.remove(file)
                print(f'Removed intermediate file: {file}')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='path', required=True, help='The path to input FASTQ or FASTQ.gz file')
    parser.add_argument('-t', metavar='threads', required=False, default=1,
                        help='The number of threads used to run the BWT construction.')
    parser.add_argument('--save_intermediate', action='store_true', help='Save intermediate files')
    args = parser.parse_args()
    save_intermediate = args.save_intermediate
    num_threads = int(args.t)
    if num_threads <= 0:
        raise ValueError(f"Number of threads must be a positive integer. You used {num_threads}.")
    main(args.i, num_threads, save_intermediate)

