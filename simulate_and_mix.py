import os
import time
import shutil
import argparse
import fastqmixer
from Bio import SeqIO


def find_genome_length(fasta_file):
    total_len = 0
    reader = SeqIO.parse(fasta_file, 'fasta')
    for contig in reader:
        total_len += len(contig.seq)
    return total_len


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-bf', '--base_fasta',
                        type=str,
                        required=True,
                        help='Path to fasta-formatted file for base genome.')
    parser.add_argument('-cf', '--contaminant_fasta',
                        type=str,
                        required=True,
                        help='Path to fasta-formatted file for base genome.')
    parser.add_argument('-d', '--depth',
                        type=float,
                        default=60,
                        help='Coverage depth desired for output genome. Defaults to 60X.')
    parser.add_argument('-t', '--tmpdir',
                        type=str,
                        default=str(time.time()).split('.')[0],
                        help='Temporary directory name.')
    parser.add_argument('-fc', '--fraction_contamination',
                        type=float,
                        required=True,
                        help='Contamination fraction. Must be between 0 and 1.')
    parser.add_argument('-o', '--output_directory',
                        type=str,
                        default=os.getcwd(),
                        help='Output directory for your FASTQ files. Defaults to current working directory.')
    args = parser.parse_args()

    if not os.path.isdir(args.tmpdir):
        os.makedirs(args.tmpdir)

    # Create FASTQs for contaminant and base genome, put them in the tmpdir for this run.
    fastqmixer.create_fastq_from_fasta(fasta_file=args.contaminant_fasta,
                                       output_fastq=os.path.join(args.tmpdir, 'contamination'),
                                       depth=args.depth)
    fastqmixer.create_fastq_from_fasta(fasta_file=args.base_fasta,
                                       output_fastq=os.path.join(args.tmpdir, 'base'),
                                       depth=args.depth)
    genome_size = find_genome_length(args.base_fasta)
    fastqmixer.mix_fastqs(genome_size=genome_size,
                          base_genome=os.path.join(args.tmpdir, 'base'),
                          contaminant_genome=os.path.join(args.tmpdir, 'contamination'),
                          depth=args.depth,
                          fraction_contamination=args.fraction_contamination)
    final_output_name = os.path.split(args.base_fasta)[-1].split('.')[0] + '_' + os.path.split(args.contaminant_fasta)[-1].split('.')[0]
    final_output_name = os.path.join(args.output_directory, final_output_name)
    cmd = 'mv sim_contam_R1.fastq.gz {output}_{fc}_R1.fastq.gz'.format(output=final_output_name,
                                                                       fc=args.fraction_contamination)
    os.system(cmd)
    cmd = 'mv sim_contam_R2.fastq.gz {output}_{fc}_R2.fastq.gz'.format(output=final_output_name,
                                                                       fc=args.fraction_contamination)
    os.system(cmd)
    shutil.rmtree(args.tmpdir)
