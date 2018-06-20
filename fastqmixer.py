import os
import shutil
import subprocess
from biotools import bbtools


def create_fastq_from_fasta(fasta_file, output_fastq, depth=20, read_length=250, insert_size=350, insert_std=10,
                            platform='MSv1'):
    # Use ART to simulate us some reads
    cmd = 'art_illumina -ss {platform} -i {input_fasta} -l {read_length} -na -p -f {depth} -m {insert_size}' \
          ' -s {insert_std} -o {output_fastq}'.format(input_fasta=fasta_file,
                                                      output_fastq=output_fastq,
                                                      depth=depth,
                                                      read_length=read_length,
                                                      insert_size=insert_size,
                                                      insert_std=insert_std,
                                                      platform=platform)
    subprocess.call(cmd, shell=True)
    # Rename the files and gzip them for space savings.
    cmd = 'mv {fastq_name}1.fq {fastq_name}_R1.fastq && gzip {fastq_name}_R1.fastq'.format(fastq_name=output_fastq)
    os.system(cmd)
    cmd = 'mv {fastq_name}2.fq {fastq_name}_R2.fastq && gzip {fastq_name}_R2.fastq'.format(fastq_name=output_fastq)
    os.system(cmd)


def mix_fastqs(genome_size, base_genome, contaminant_genome, fraction_contamination=0.2, depth=20, tmpdir='tmp'):
    # Make tmpdir if it doesn't already exist.
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    # Subsample our base genome to desired coverage level.
    # BBTools think we want to sample everything if num bases is set to 0, so don't do if that is the case.
    if 1.0 - fraction_contamination != 0.0:
        bbtools.subsample_reads(forward_in=base_genome + '_R1.fastq.gz',
                                forward_out=os.path.join(tmpdir, base_genome + '_R1.fastq.gz'),
                                reverse_in=base_genome + '_R2.fastq.gz',
                                reverse_out=os.path.join(tmpdir, base_genome + '_R2.fastq.gz'),
                                num_bases=str(genome_size * (1.0 - fraction_contamination) * depth))
    # Repeat process for contaminant genome
    if fraction_contamination != 0.0:
        bbtools.subsample_reads(forward_in=contaminant_genome + '_R1.fastq.gz',
                                forward_out=os.path.join(tmpdir, contaminant_genome + '_R1.fastq.gz'),
                                reverse_in=contaminant_genome + '_R2.fastq.gz',
                                reverse_out=os.path.join(tmpdir, contaminant_genome + '_R2.fastq.gz'),
                                num_bases=str(genome_size * fraction_contamination * depth))
    # Mix together the two sets of subsampled reads.
    cmd = 'cat {base_genome_forward} {contaminant_genome_forward} > ' \
          '{output_filename}'.format(base_genome_forward=os.path.join(tmpdir, base_genome + '_R1.fastq.gz'),
                                     contaminant_genome_forward=os.path.join(tmpdir, contaminant_genome + '_R1.fastq.gz'),
                                     output_filename='sim_contam_R1.fastq.gz')
    os.system(cmd)
    cmd = 'cat {base_genome_forward} {contaminant_genome_forward} > ' \
          '{output_filename}'.format(base_genome_forward=os.path.join(tmpdir, base_genome + '_R2.fastq.gz'),
                                     contaminant_genome_forward=os.path.join(tmpdir, contaminant_genome + '_R2.fastq.gz'),
                                     output_filename='sim_contam_R2.fastq.gz')
    os.system(cmd)
    shutil.rmtree(tmpdir)


