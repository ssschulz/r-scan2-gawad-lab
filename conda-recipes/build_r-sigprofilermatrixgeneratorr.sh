#!/bin/bash

conda-build -c conda-forge -c bioconda -c jluquette r-sigprofilermatrixgeneratorr

anaconda upload \
    /home/ljl11/miniconda3_for_building/conda-bld/linux-64/r-sigprofilermatrixgeneratorr-1.0-py38r41h3fd9d12_0.tar.bz2
