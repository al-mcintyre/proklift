# proklift
A liftover program to convert genomic indices between strains or species based on multiple genome alignment from Mauve (in xmfa format from gui - primary output, see http://darlinglab.org/mauve/user-guide/files.html)
```
usage: proklift.py [-h] -x ALIGNMENT -b BED [-o OUTPUT] [-v] 
```
arguments:
```
  -h, --help            show this help message and exit
  -x ALIGNMENT, --alignment ALIGNMENT
                        alignment file (xmfa format) from Mauve GUI for two
                        species
  -b BED, --bed BED     bed file with positions to exchange - if no strand
                        column provided, chooses + by default for output
  -o OUTPUT, --output OUTPUT
                        output file name
  -v, --version         print version
```
output format:
bed file with 
second genome "2",start,end,whether strand matches input (0/1),".",strand, 
followed by first genome indices: "1",start,end  
