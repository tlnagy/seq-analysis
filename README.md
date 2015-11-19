# seq_analysis

ET0H's analysis pipeline for the ubiquitin sequencing data. 

## Usage

Ipython notebook is recommended for interacting with the data. Assuming
you are in the parent directory of `seq-analysis`, use as follows:

```python
import pandas as pd
aa_weighted = pd.read_csv("seq-analysis/data/aa_weighted.csv", index_col=0)
# convert to Multiindex
aa_weighted.set_index(["days", "amino acids", "positions"], append=True, inplace=True)
# getting only the Et0H weighted mean slope data and have the days in 
# separate columns
aa_weighted.loc["Et0H", "weighted mean slope"].unstack("days")
```

which looks like:

## Setup

### Code

We're using Anaconda to manage all packages. First install Miniconda from
<http://conda.pydata.org/miniconda.html>. Either Python 2 or Python 3 will
work, because we will be creating a virtual environment. It should add the
`conda` executable to your path, if `which conda` give you nothing then
make sure to add `conda`. Then run

```
conda create -n seqanalysis -y python=3 seaborn ipython-notebook biopython scipy scikit-learn
source activate seqanalysis
ipython notebook
```

### Data

The raw processed data is called `lanes_new_combined.fastq.h5` and it's
located in the `et0h/data` folder on the derisilab105 server.

## Structure

There are two main files: `seq_analysis.py` and `demultiplex_fastq.py`.
The latter takes a raw fastq file and demultiplexes it. The former
actually goes through and processes the data. Most of the interesting
stuff is there.
