Mining association rules in KEGG orthologs of human gut bacteria
==================================================================
The project was done in Tel-Aviv University. <br />
Personal instructor: Dr. Efrat Muller <br />
Head of Laboratory: Prof. Elhanan Borenstein <br />
School of Computer Science and Faculty of Medicine, Computational Metagenomics Research Lab 
-------------------------------------------------------------------

This repository includes the code for the R script and R notebooks I used during my research under Dr. Muller. 

To those interested, this is the paper itself:
> [Mining association rules in KEGG orthologs of human gut bacteria](Alpha-Project/Final_Project_To_Upload.pdf)

The goal of our research was two-pronged - provide a proof of concept to the effectiveness of using association rule mining to research the microbiome, and try to find interesting biological results with said association rules. 

We succesfully used the `apriori` algorithm with `hyperedge sets` to mine association rules, and the following is a taste of our results: 

![image](https://github.com/user-attachments/assets/ff6ad37b-4784-4000-95d4-b91c1522848d)

We used different techniques to narrow down the rules to some non-trivial ones, and found a few with potential to be important:

> 1. {K03091,K06409,K06960}: According to KEGG database, K03091 and K06409 are both
genes tied to sporulation in certain bacteria, while K06960 is an undocumented protein. All
three of them appear in 500 different bacterial species, but K06960 appears almost
exclusively in bacteria in which K03091 and K06409 appear. That means that there’s a
chance it’s tied to sporulation, or another function performed exclusively by bacteria that
can perform sporulation. Additionally, K06960 has been found to be predictive of lynch
syndrome, which is the most common cause of hereditary colorectal cancer, meaning that
in patients in which it appears there is a high likelihood of developing colorectal cancer
(Yan et al., 2020). There could also be a connection between that and the other two genes in
the set, but more extensive research will need to be done.
> 2. {K18888,K18887,K16785,K07166,K09157}: According to KEGG database, the first three
KO’s - {K18888,K18887,K16785} – all appear within a single metabolic pathway, pathway
2010 – which is a pathway of ABC transporters, which are proteins that are responsible for
translocating solutions using the energy from ATP hydrolysis (Jones & George, 2004).
However, the second two KO’s – {K07166,K09157} - aren’t uncharacterized in the
literature, and we suggest that they might also be related to ABC transporters.
