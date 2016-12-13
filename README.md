# *in vivo* RBP motif prediction by target site accessibility #



## Input file##

## Pre-requirement ##

RNAPLFOLD (version 1.7.2)
This package requires accessibility calculation in advance, which can be computed by RNAplfold. For detailed information of running RNAplfold, Please refer to http://www.tbi.univie.ac.at/~ivo/RNA/RNAplfold.html 

For running our scripts, it should be sufficient to do (assuming the length of the binding site is 8 nt) :

RNAplfold -W 80 -L 40 -u 8  < seq_file_name  > /RNAplfold_direct /W80L40u8_whole.txt 

(The setting of parameters W, L, u can be changed upon different situation, but please make sure the output file name is also changed to reflect the parameters’ setting.)


## Related Publications ##

X. Li, G. Quon, H.D. Lipshitz, Q.D. Morris, **Predicting in vivo binding sites of RNA-binding proteins using mRNA secondary structure**, *RNA* 16.6 (2010): 1096–1107. [[Pubmed]](https://www.ncbi.nlm.nih.gov/pubmed/20418358)
