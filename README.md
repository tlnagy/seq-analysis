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

group_name = "Et0H"
_, mapped_barcode_data = seq.map_dataset("seq-analysis/lanes_new_combined.fastq.h5", group_name)
processed = seq.process_data(mapped_barcode_data)
regressed = seq.regress(processed, group_name, "seq-analysis/truseq_primers.csv")
filtered = seq.groupby_filter(regressed)
filtered
```

## Structure

There are two main files: `seq_analysis.py` and `demultiplex_fastq.py`.
The latter takes a raw fastq file, demultiplexes it, and produces
something like `et0h_barcodes_to_count.csv`. The former actually goes
through and processes the data. Most of the interesting stuff is there.
