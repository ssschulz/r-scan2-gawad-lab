R package with supporting routines for the SCAN2 pipeline.  These routines
may not be particularly useful when separate from the main pipeline.

See https://github.com/parklab/SCAN2 for the full pipeline.


# Major caveat
This package requires command line access to a Linux system. There are 3
reasons:

1. Custom tabix parsing to read in large files. This uses a command line
   `tabix` executable piped into the \*nix `cut` command. Sometimes tabix is
   replaced by a simple `gunzip -c` when no genomic region is given.
2. The bedtoolsr package: this is a wrapper around bedtools executables.
3. The SigProfilerMatrixGeneratorR package: another wrapper (via the R
   R<->python reticulate interface) to a python library.


# Extra installation steps
After installing r-scan2, the user must download and install a reference
genome for SigProfilerMatrixGenerator.  See the SCAN2 pipeline README
for detailed instructions.
