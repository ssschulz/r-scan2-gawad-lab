#!/bin/bash

conda-build -c conda-forge -c bioconda -c jluquette r-scan2

anaconda upload  /home/ljl11/miniconda3_for_building/conda-bld/linux-64/r-scan2-1.0-py38r41h3fd9d12_4.tar.bz2
