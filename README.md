#ConCat v1.0

###BioPython and Blast command line is a prerequisite for Codon-Align package

##Requirements:

Python 2.7 and above

In the firt step Codon-Align extracts all trasncript sequences for a selected species group and stores it in Sequences folder

``` 
python fetch.py -cds aa.txt -ortho
```
In the second step Condon-Align uses CCDS exon annotations to perform exon wide search accross all the species gene trasncript sequences via tblastn search. Individual exon sequences are aligned using muscle codon alignment method and finally concatenated together to build final gene codon sequence alignemnt.  

``` 
python python consensus.py -i aa.txt -o OutputFolder
```



