# logan_unsticker
"Sticky sequences" are sequences that appear in common genomes on the SRA (human, mouse, yeast etc...). If an input protein database contains these sequences (or subsequences), the output of the DIAMOND search of Logan will be highly populated with false positive contigs that have been captured by these sticky sequences. 

logan_unsticker.sh is a short script designed to identify "sticky sequences" in a given protein database against a common genome of choice. For example, if generating a papillomavirus database from NCBI nucleotide sequences, this script can identify which (and where) sequences contain human host DNA from papillomavirus integration sites.

See Issues tab for recommended genomes to screen against - genomes to be included will be continuously updated. False positives can be commented here for inclusion.

![image](https://github.com/syueqiao/logan_genome_QC/assets/105825554/5fb5cf9a-1d70-4d29-8dd5-adac3b96ef2b)

### Inputs
The script as it is set up currently accepts 2 inputs:
1. A FTP url of the genome to be scanned, either in 2bit, or gzipped fasta format. This can be easily manipulated to accept other input formats.
2. A fasta file containing protein sequences to generate a diamond database with. Script can also be modified to directly accept formatted diamond .dmnd databases.

To toggle between 2bit/fasta format, use -b flag to denote 2bit format. 
twoBitToFa utility can be found: https://hgdownload.soe.ucsc.edu/admin/exe/

### Outputs
The script outputs a diamond .pro file, which can be used to remove whole or partial sequences depending on desired protein DB and context. 

### Example Usage
With  2bit file format:
```
./logan_unsticker.sh -b "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit" PV_test_DB.fasta
```
With a gzipped fasta file:
```
./logan_unsticker.sh "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" PV_test_DB.fasta
```
### Example Output
in sacCer3.pro
```
NC_000913.3_sliding:839201-840200	411	73	1000	-	papilloma.Late_protein_L1.Human_papillomavirus:MH777234	263	376	518	27.9	6.90e-04	6M5I32M2I25M3D9M1I18M6D15M	IKNGKDVAQDGSSSLPYTPEHAFTLWSQYQATDDISVGAGARYIGSMHKGSDGAVGTPAFTEGYWVADAKLGYRVNRNLDFQLNVYNLFDTDYVASINKSGYRYHPGEPRTFL	ACATCATACAATTTCTCCAAAAAGTGGGGCCTGCGCCCCACATCTGAATCAGAAATGCATATTGGCTGTGAGCAAGAAGGTTCTTGGCTCGCCCGGGTGATAACGGTAGCCGCTCTTGTTGATTGAGGCGACGTAATCGGTATCAAACAGGTTGTAAACGTTTAGCTGGAAGTCGAGATTGCGATTAACTCGATACCCCAGTTTGGCATCGGCGACCCAGTAACCTTCGGTAAACGCTGGCGTTCCCACCGCGCCGTCTGAACCTTTATGCATACTGCCGATATAGCGTGCGCCCGCGCCAACAGAGATATCGTCGGTTGCCTGATATTGGCTCCATAAGGTGAAGGCGTGCTCCGGGGTATACGGCAGCGATGAGGAACCATCCTGGGCAACATCTTTGCCGTTTTTGATGGTTGCTTTTTGCTGGGTATAGCCGCCAATCACCTGCCACGCGGGAGTGATATTCCCGGCCACGGATATCTCATAGCCTTCGACGCGTTTCTTACCGTATTGCGAGTAAGTTCCGTCATCATTTTGCTCAACTTCATTTTCGATATCAGTGCGGAACAGCGCGGCGGTGAGCAACAGACGTTTATCCAGAACCTGCCATTTGGTGCCAATCTCGCTGGTGTTGGCTTTTTGCGGTTTAAAATCGGTGCGGTTGGCACTGTTACCGCTGCCAGACTGCGCAAGGGCGAAGTTGTTGCCGCCCGGAGGCTGCTGGGAAACGGCATAGTTAATATAGACATTGCCGTTTTCCGTCAGGTGATACAGCGCCCCGGCTTTCCAGTTCATCAGATTGCCCGACTTGGCGGTGTCGACGGTGGTGACCGGAGAACCTTTTGCCACACCAGTTGGGCAGGTGATGGCACCGCGTCCGCTGCCGCCGCAGGCGGTGGCACTGTCATATTCAGTATGATAATTATCCAGACGGATCCCGCCGTTCAGCTCAAAATCACGGGTGATTTGCAGCGTATCAAAGGCGTAAATTGCGAAGGTA	*
NC_000913.3_sliding:3988801-3989800	914	666	1000	-	papilloma.Late_protein_L2.Papillomavirus_panthera7600:BK066462	414	495	525	32.2	9.13e-04	8M1I9M2I19M2I5M4D37M	QLPQLADQLAALGESDLLFALSQHAVAFAQSQLHQQDRKWPRLPDYFAIGRTTALALHTVSGQKILYPQDREISEVLLQLPEL	TCGGCTTTTTGGCTCTCCTGGGCTTTTTGCAATGCCGTCAGTTGGTTAGCCAGGGCATCGCTGGTGGCGGTCTGATTGACGGCCTGTTGTTTACCCCAGCCATACAAACCGATGCCCGCCGCCAGAGCAATAGCGATAGCCACCGCGCTGAGAATCAATGCGGTATTGTTCTTACTCTTTTTTTCTGTTGCGACAGGTTGTGACGTGGTGTCCACGGCCTCCCTGGTCTCTTCAACCACGGCGGAGGTTTTTTCTTGTTCCGTCATTATGGCTTCCTGTTATGAGAGTTATTGTAATGCCCGTAAAAGCGCATCGTTGTCAGCGTTATCGGCGACCTTAATGTCTTGCCAGCCCAGTTCCCGGGCGAGTTTCGCCAAACGCTCACTGACGACCAATAGTCGACAGTGTAGTAACCAGTGCTCACGATACCATTGTGGGATCAGCGACCAGAGTTGCTGCAACATTTCACCGCTGGTAACAACGACCATCGTCACCTCGCGGGCTTGCCAGCGCATCGCTTCTTCTGCACCATCGTAATGGATTGCGCATCGTTGATAACATTCACAAAAAGTGACCTCAGCACCGCGCGCCGTCAGGGTATCCCCAATTAGCTCACGACCACCATTGCCACGTAATATCAGCGCACGTTTGCCCGCAATATTTTGTAATTCAGGTAATTGTAGCAAGACTTCGCTGATTTCCCGATCCTGCGGGTAGAGAATCTTCTGTCCACTTACGGTATGTAGTGCCAGTGCGGTGGTGCGTCCAATGGCGAAATAATCAGGTAGTCGGGGCCATTTACGATCTTGCTGATGCAGCTGTGATTGGGCAAAAGCAACCGCGTGTTGCGAGAGGGCAAACAACAGATCGCTCTCCCCCAGCGCTGCCAGTTGATCAGCAAGTTGCGGTAATTGTTGACCCGGAGAAAACTCAATCAGCGGAAAATGCCAGGCCACCTGCCCCAGTGTGCGCAGACGGCTCACTAACTCTTCTCCAGCGG	*
NC_000913.3_sliding:4339201-4340200	319	5	1000	-	papilloma.Late_protein_L1.Bandicoot_papillomatosis_carcinomatosis_virus_type_2:EU277647	237	339	506	28.0	5.45e-05	9M2I23M1D20M2I20M1D29M	DLLEKTSDRLHFDEAWYGYARFNPIYADHYAMRGEPGDHNGPTVFATHSTHKLLNALSQASYIHVREGRGAINFSRFNQAYMMHATTSPLYAICASNDVAVSMMD	TGCCGTCCATCATCGACACCGCCACGTCGTTGGATGCGCAGATGGCATACAGCGGGGAGGTGGTGGCATGCATCATGTAGGCCTGGTTGAAGCGGGAGAAGTTAATCGCCCCACGACCTTCACGTACATGAATATAAGAAGCCTGTGACAGCGCATTCAGCAGTTTGTGGGTGGAGTGGGTGGCGAAAACGGTAGGACCGTTGTGATCGCCAGGTTCGCCGCGCATGGCATAGTGATCGGCATAGATCGGGTTGAAACGTGCATAGCCGTACCAGGCTTCGTCAAAGTGCAGACGATCGGAGGTTTTTTCCAGCAGATCCTGCGCTTCTTTAGCGTTATAACACACGCCGTCATAGGTGCAGTTGGTCACCACGCAGTAAGACGGTTTTTGCCCGGCTTTGTCTTTGGTCAGCGGGCTTTCACTGATTTTCTTCTGCAAGGTTTCAGGTTGCATTTCCTGCGGATAGATTGGCCCGATAATGCCGTAGCGGTTGCGGCTTGGCACCATATAGACCGGTTTCGCGCCTGTCAGCATCAAACCTTGTTCGATGGATTTATGGCAGTTACGGTCAACGACCACGACATCGTTATCGGTCATGCAAGCCTGCATGATGGTGCGGTTAGAGCCGGAAGTACCGACGACTACCGACCAGGAGCGATCGGCACCAAATACGCGTGCGGCATATTTTTCGCTTTCGCCAAATGCGCCAGTATGGTCAAGCAAAGAACCGAGGGAAGTTCGTTCGATGCCCATGTCGGTGCGGAACAGATTTTCACCATAGTAGTCATGGTAGAAACGTCCGGCGGGTGTTTTGGTAAAACCAACGCCGCCCTGGTGGCCTGGCGCTGCCCAGGAATATTCATGGATGTCACTATATTTCATCAGCGCGCTGAACAGTGGCGGCAACAGCTGCTGGCGGTAGCGGGTCATCGCGGCAACGGCGCGTCCGGCGATAAAGTCGGCGGTATCTTCCAGAATCCAGGCGAATTCATCGACAAG	*
```
