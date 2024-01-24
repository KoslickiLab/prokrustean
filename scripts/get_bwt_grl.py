import os
import glob
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
def main(input):
    # input file fastq name ex. "./downloads/x_x_x.fastq.gz"
    in_path_for_fastq_file = input
    
    if str.endswith(input, '.fastq.gz'):
        extension='.fastq.gz'
    elif str.endswith(input, '.fastq'):
        extension='.fastq'
    else:
        raise Exception('non-fastq file (fastq.gz or .fastq)')

    # input file becomes concatenated ex. "./downloads/x_x_x.txt"
    in_path_for_concatenated_string = in_path_for_fastq_file.replace(extension, ".txt")
    # grlbwt processed. they don't like underscore ex. "./downloads/xxx.txt"
    out_path_grlbwt_rl_bwt_file = in_path_for_fastq_file.replace(extension, "").replace("_", "")
    # grlbwt processed to text. ex. "./downloads/x_x_x.bwt"
    out_txt_path_grlbwt_txt_file = in_path_for_fastq_file.replace(extension, ".bwt")
    print(f'input: {input}')
    print(f'sequences will be concatenated : {in_path_for_concatenated_string}')
    print(f'grlbwt will be processed : {out_path_grlbwt_rl_bwt_file}')
    print(f'grlbwt output will be converted to txt : {out_txt_path_grlbwt_txt_file}')

    def process_fastq_gz(file_path, output_file_path, extension):
        concatenated_sequence = ""
        if extension=='.fastq.gz':
            with gzip.open(file_path, "rt", encoding="utf-8") as file:  # 'rt' stands for read text
                lines = file.readlines()
                for i, line in enumerate(lines):
                    if (i % 4 == 1):  # sequences are located every 4 lines starting at line 1 (0-indexed)
                        concatenated_sequence += line.strip() + '$'
        elif extension=='.fastq':
            with open(file_path, "rt") as file:  # 'rt' stands for read text
                lines = file.readlines()
                for i, line in enumerate(lines):
                    if (i % 4 == 1):  # sequences are located every 4 lines starting at line 1 (0-indexed)
                        concatenated_sequence += line.strip() + '$'
        else:
            raise Exception('non-fastq extension')
        with open(output_file_path, 'w') as file:
            file.write(concatenated_sequence)

    def run_grlbwt(file_path, out_file_path, out_txt_file_path):
        os.system(f'../../grlbwt/build/grlbwt-cli {file_path} -o {out_file_path}')
        os.system(f'../../grlbwt/build/grl2plain {out_file_path}.rl_bwt {out_txt_file_path}')

    # preprocess 
    if os.path.isfile(in_path_for_concatenated_string):
        print('already concatenated')
    else:
        process_fastq_gz(in_path_for_fastq_file, in_path_for_concatenated_string, extension)

    # grlbwt and conversion to txt 
    run_grlbwt(in_path_for_concatenated_string, out_path_grlbwt_rl_bwt_file, out_txt_path_grlbwt_txt_file)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', metavar='path', required=True, help='the path to input')
    args = parser.parse_args()
    main(args.i)

