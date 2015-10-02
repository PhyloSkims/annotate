
Due to the long compilation time of ncbi blast sources (~1h),
the default build actually consists in extracting pre-compiled
binaries.

To revert to full compilation from sources, please edit
the Makefile and change :

DIRS = binaries

to

DIRS = src


