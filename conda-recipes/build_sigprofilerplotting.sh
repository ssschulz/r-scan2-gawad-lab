#!/bin/bash

source /home/ljl11/miniconda3_for_building/bin/activate 

grayskull pypi --strict-conda-forge sigprofilerplotting

# Add numpy < 1.23.0 requirement (currently numpy is 1.23.1, causes scipy to fail)
awk '{ if ($0 ~ /- python/) { print $0; print "    - numpy <1.23.0"; } }' sigprofilerplotting/meta.yaml > sigprofilerplotting/tmpmeta.yaml
mv sigprofilerplotting/tmpmeta.yaml sigprofilerplotting/meta.yaml 


conda-build -c conda-forge -c bioconda sigprofilerplotting


anaconda upload     /home/ljl11/miniconda3_for_building/conda-bld/noarch/sigprofilerplotting-1.2.2-py_0.tar.bz2
