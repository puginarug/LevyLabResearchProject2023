here you can find the scripts and commandline promps used for the local analyses in this work.

- Colabfold - for the folding of toxin and immunity proteins.
- Foldseed - for localy comparing the structure of proteins.
- Code used for scooring the interaction between toxin and immunity:
-- analysepairs.py script used to calculate the mean intermolecular PAE score (MIP) for each pair (adapted from:   https://github.com/alexlevylab/T6E_discovery/blob/main/1_scoring_of_match_upload.py)
-- a folder containig the logistic regression model used to determine weather there is an interaction based on the MIP.
-- a link to the pdockq.py script used to calculate the pDockQ score: https://gitlab.com/ElofssonLab/FoldDock/-/blob/main/src/pdockq.py
