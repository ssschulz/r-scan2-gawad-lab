#!/bin/bash

source /home/ljl11/miniconda3_for_building/bin/activate 

echo "================================================================================="
echo "IMPORTANT"
echo "IMPORTANT: if you've just built sigprofilerplotting you must avoid using the locally built package"
echo "IMPORTANT: One way is to `anaconda upload` it and then delete or move the built package at "
echo "IMPORTANT: /home/ljl11/miniconda3_for_building/conda-bld/noarch/sigprofilerplotting-1.2.2-py_0.tar.bz2"
echo "IMPORTANT: The package file cannot just be renamed in that directory or else conda will automatically"
echo "IMPORTANT: find it (the package name is stored in the tarball)"
echo "IMPORTANT"
echo "================================================================================="
if [ -f "/home/ljl11/miniconda3_for_building/conda-bld/noarch/sigprofilerplotting-1.2.2-py_0.tar.bz2" ]; then
    echo "see above warning; please remove this file"
    exit 1
fi


# SigProfilerMatrixGenerator is a package on PyPI. conda provides an
# automatic way to construct conda packages with `conda skeleton pypi`.
echo "================================================================================="
echo "conda skeleton pypi"
echo "================================================================================="
conda skeleton pypi --python-version 3.8 --version 1.2.9 SigProfilerMatrixGenerator

# This is fun: the source tarball for sigprofilermatrixgenerator has a case-sensitive URL.
# Yep. Really. This fails (always, not intermittent) with a 404 error:
#     wget https://pypi.io/packages/source/s/sigprofilermatrixgenerator/sigprofilermatrixgenerator-1.2.9.tar.gz
# and this works:
#     wget https://pypi.io/packages/source/s/sigprofilermatrixgenerator/SigProfilerMatrixGenerator-1.2.9.tar.gz
# The above issue is fixed by using SigProfilerMatrixGenerator in place of sigprofilermatrixgenerator
# in the conda skeleton pypi command above.

echo "Removing auto-generated test imports that DO NOT WORK even on pip-installed SigProfilerMatrixGenerator:"
echo "    import SigProfilerMatrixGenerator.references.chromosomes.tsb"
echo "    import SigProfilerMatrixGenerator.references.matrix"
echo "    import SigProfilerMatrixGenerator.references.vcf_files.BRCA_example.SNV"

sed -i -e 's/\(- SigProfilerMatrixGenerator.references.chromosomes.tsb\)$/#\1/' -e 's/\(- SigProfilerMatrixGenerator.references.matrix\)$/#\1/' -e 's/\(- SigProfilerMatrixGenerator.references.vcf_files.BRCA_example.SNV\)$/#\1/' sigprofilermatrixgenerator/meta.yaml

# -c jluquette makes our local sigprofilerplotting package available
# IMPORTANT: be sure tha
echo "================================================================================="
echo "conda-build"
echo "================================================================================="
conda-build -c conda-forge -c bioconda -c jluquette sigprofilermatrixgenerator


echo "================================================================================="
echo "anaconda upload"
echo "================================================================================="
anaconda upload /home/ljl11/miniconda3_for_building/conda-bld/linux-64/sigprofilermatrixgenerator-1.2.9-py38_0.tar.bz2
