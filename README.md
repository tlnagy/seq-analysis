# seq_analysis

ET0H's analysis pipeline for the ubiquitin sequencing data. 

## Usage

Ipython notebook is recommended for interacting with the data. Run as
follows:

```python
%load_ext autoreload
%autoreload 2

import sys
sys.path.append('seq-analysis/')
import seq_analysis as seq

slopes = seq.process_data("seq-analysis/lanes_new_combined.fastq.h5")
slopes
```

## Structure

There are two main files: `seq_analysis.py` and `demultiplex_fastq.py`.
The latter takes a raw fastq file, demultiplexes it, and produces
something like `et0h_barcodes_to_count.csv`. The former actually goes
through and processes the data. Most of the interesting stuff is there.
