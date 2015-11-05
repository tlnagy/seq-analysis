# seq_analysis

ET0H's analysis pipeline for the ubiquitin sequencing data. 

## Usage

Ipython notebook is recommended for interacting with the data. Run as
follows:

```
%load_ext autoreload
%autoreload 2

import seq_analysis as seq
df = seq.load_dataset("et0h_barcodes_to_count.csv")
df
```
