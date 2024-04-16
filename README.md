# genome_info_

[![Build Status](https://travis-ci.org/kylessmith/genome_info.svg?branch=master)](https://travis-ci.org/kylessmith/genome_info) [![PyPI version](https://badge.fury.io/py/genome_info.svg)](https://badge.fury.io/py/genome_info)
[![Coffee](https://img.shields.io/badge/-buy_me_a%C2%A0coffee-gray?logo=buy-me-a-coffee&color=ff69b4)](https://www.buymeacoffee.com/kylessmith)

A package to work with genomic intervals.


## Install

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  
```
    https://continuum.io/downloads
```

Dependencies can be installed by:

```
    pip install -r requirements.txt
```

PyPI install, presuming you have all its requirements installed:
```
    pip install genome_info(not yet!!)
```

## Usage

```python
import genome_info_
genome = genome_info.GenomeInfo("hg38")
genome.chrom_sizes

genome.n_bases

genome.tss()

genome.tss(upstream=1000, downstream=1000)

```

