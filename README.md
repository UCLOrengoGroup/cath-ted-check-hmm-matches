# Check TED/HMM matches

Check domain boundaries between HMM matches and TED consensus predictions.

The incoming data is a subset of all HMM/TED domains corresponding to UniProt accessions from human.

Currently the processing is done in a Jupyter notebook which:
* filters the incoming data based on minimum domain length 
* creates bins according to overlap and evalue
* randomly selects 10 entries per bin  

Ideally this would also include an easy way of visualising the domain boundaries and recording whether the domain boundaries are good/bad.
