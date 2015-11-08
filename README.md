# seq_analysis

ET0H's analysis pipeline for the ubiquitin sequencing data. 

## Usage

Ipython notebook is recommended for interacting with the data. Run as
follows:

```python
%load_ext autoreload
%autoreload 2

import seq_analysis as seq
_, mapped_barcode_data = seq.map_dataset("et0h_barcodes_to_count.csv")
processed_barcodes = seq.process_data(mapped_barcode_data)
df = seq.regress(processed_barcodes)
df
```

## Structure

There are two main files: `seq_analysis.py` and `demultiplex_fastq.py`.
The latter takes a raw fastq file, demultiplexes it, and produces
something like `et0h_barcodes_to_count.csv`. The former actually goes
through and processes the data. Most of the interesting stuff is there.
