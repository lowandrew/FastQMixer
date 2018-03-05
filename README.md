# FastQMixer

Given some FASTA files, will simulate FASTQs from them using ART, and then mix them together in the desired proportions.
Should work on any UNIX-based system.

# Dependencies

- The ART executable (art\_illumina) must be on your $PATH.

# Python Packages

- biopython >= 1.70
- OLCTools >= 0.3.45

To install both, use pip: `pip install biopython olctools`

# Usage

Clone this repository - you'll want to use the `simulate_and_mix.py` script.

```
usage: simulate_and_mix.py [-h] -bf BASE_FASTA -cf CONTAMINANT_FASTA
                           [-d DEPTH] [-t TMPDIR] -fc FRACTION_CONTAMINATION

optional arguments:
  -h, --help            show this help message and exit
  -bf BASE_FASTA, --base_fasta BASE_FASTA
                        Path to fasta-formatted file for base genome.
  -cf CONTAMINANT_FASTA, --contaminant_fasta CONTAMINANT_FASTA
                        Path to fasta-formatted file for base genome.
  -d DEPTH, --depth DEPTH
                        Coverage depth desired for output genome.
  -t TMPDIR, --tmpdir TMPDIR
                        Temporary directory name.
  -fc FRACTION_CONTAMINATION, --fraction_contamination FRACTION_CONTAMINATION
                        Contamination fraction. Must be between 0 and 1.

```


