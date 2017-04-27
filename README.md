# An RBP motif prediction algorithm based on target site accessibility #

This motif discovery algorithm selects the final motif with most distinguished feature (total accessibility or site number) in bound (positive) and unbound (negative) transcripts. 

## Pre-requirements ##

* [Python 2.7](https://www.python.org/downloads/)

* [RNAplfold (version 2.1.8)](https://www.tbi.univie.ac.at/RNA/index.html)

    Sequence accessibility of transcripts is predicted by the RNAplfold tool in the ViennaRNA package.

## Input files ##

* **pos_file**: A file that contains gene names in the positive set, one gene per line.

* **neg_file**: A file that contains gene names in the negative set, one gene per line.

* **seq_file**: A fasta file containing sequences of genes in positive and negative sets (in the __RNA alphabet__).

* **RNAplfold_direct**: A folder containing the `<ID>_lunp` files of RNAplfold output of genes in positive and negative sets. Make sure gene names in this file are consistent with those in pos_file and neg_file.

    To get result from RNAplfold, assuming the length of the binding site is at most 10 nt:

    ```
    RNAplfold -W 80 -L 40 -u 10  < seq_file_name
    ```

    The setting of parameters W, L, u can be changed upon different situation.

    For detailed information of running RNAplfold, please refer to [https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html](https://www.tbi.univie.ac.at/RNA/RNAplfold.1.html).


## Parameters ##

To run the code, you will need to specify the following parameters by modifying the following lines of `run2_new_RNAplfold_format.py`:

```
model = 'access_seq'
pos_file_name = ''   # name of pos_file
neg_file_name = ''   # name of neg_file 
seq_file_name = ''   # name of seq_file
RNAplfold_direct = ''  # directory name of the RNAplfold results

final_out_name = ''   # name of the final output (only the best motif)
detailed_final_out_name = ''   #name of the detailed output (motifs from all five seed 6mers)
```

## Motif prediction ##

To run the code, one need to run the following function in the `run2_new_RNAplfold_format.py`, 

```
Motif_discovery(final_out_name, detailed_final_out_name, model, pos_file_name, neg_file_name, seq_file_name, RNAplfold_direct)   
```

## Related Publication ##

X. Li, G. Quon, H.D. Lipshitz, Q.D. Morris, **Predicting in vivo binding sites of RNA-binding proteins using mRNA secondary structure**, *RNA* 16.6 (2010): 1096–1107. [[Pubmed]](https://www.ncbi.nlm.nih.gov/pubmed/20418358) [[Supplementary materials]](http://rnajournal.cshlp.org/content/16/6/1096/suppl/DC1)
