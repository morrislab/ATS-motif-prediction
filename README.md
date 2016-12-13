# An RBP motif prediction algorithm based on target site accessibility #

This motif discovery algorithm selects the final motif with most distinguished feature (total accessibility or site number) in bound (positive) and unbound (negative) transcripts. 

## Pre-requirement ##

* [RNAplfold (version 1.7.2)](https://www.tbi.univie.ac.at/RNA/index.html)

    Sequence accessibility of transcripts can be predicted by RNAplfold in the ViennaRNA package.
    
    For example, assuming the length of the binding site is 8 nt:

    ```
    RNAplfold -W 80 -L 40 -u 8  < seq_file_name  > /RNAplfold_direct /W80L40u8_whole.txt
    ```

    The setting of parameters W, L, u can be changed upon different situation, but please make sure the output file name is also changed to reflect the parameters’ setting.
    
    For detailed information of running RNAplfold, Please refer to [https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html](https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html).

## Input files ##

* pos_file: A file that contains gene names in the positive set, one gene per line.

* neg_file: A file that contains gene names in the negative set, one gene per line. 

* seq_file: A fasta file containing sequences of genes in positive and negative sets.

* RNAplfold_result: Result of RNAplfold output of genes in positive and negative sets. Make sure gene names in this file are consistent with those in pos_file and neg_file.

## parameters ##

To run the code, you will need to specify the following parameters by modifying the following lines of `run2_new_RNAplfold_format.py`:

```
model = 'access_seq'
pos_file_name = ''   # name of pos_file
neg_file_name = ''   # name of neg_file 
seq_file_name = ''   # name of seq_file
RNAplfold_direct = ''  # directory name of the RNAplfold results

final_out_name = ''   # name of the final output
detailed_final_out_name = ''   #name of the final output (the detailed version)
```

## Related Publications ##

X. Li, G. Quon, H.D. Lipshitz, Q.D. Morris, **Predicting in vivo binding sites of RNA-binding proteins using mRNA secondary structure**, *RNA* 16.6 (2010): 1096–1107. [[Pubmed]](https://www.ncbi.nlm.nih.gov/pubmed/20418358)
