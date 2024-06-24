# logan_genome_QC
Short script to run quality control of a given protein DB against a input genome of choice. For example, if generating a papillomavirus database, this script can be used to identify sequences that might return false positives, such as those with integration sites that may contain host or misassembled/misannotated off-target sequences. 

TBD PIPELINE IMAGE

### Inputs
The script as it is set up currently accepts 2 inputs:
1. A FTP url of the genome to be scanned, either in 2bit, or gzipped fasta format. This can be easily manipulated to accept other input formats.
2. A fasta file containing protein sequences to generate a diamond database with.

To toggle between 2bit/fasta format, use -b flag to denote 2bit format. 

### Outputs
The script outputs a diamond .pro file, which can be used to remove whole or partial sequences depending on desired protein DB and context. 

### Example Usage
With  2bit file format:
```
./logan_genome_QC.sh -b "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.2bit" PV_test_DB.fasta
```
With a gzipped fasta file:
```
./logan_genome_QC.sh "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" PV_test_DB.fasta
```
### Example Output
in sacCer3.pro
```
NC_000913.3_sliding:839201-840200	411	73	1000	-	papilloma.Late_protein_L1.Human_papillomavirus:MH777234	263	376	518	27.9	6.90e-04	6M5I32M2I25M3D9M1I18M6D15M	IKNGKDVAQDGSSSLPYTPEHAFTLWSQYQATDDISVGAGARYIGSMHKGSDGAVGTPAFTEGYWVADAKLGYRVNRNLDFQLNVYNLFDTDYVASINKSGYRYHPGEPRTFL	ACATCATACAATTTCTCCAAAAAGTGGGGCCTGCGCCCCACATCTGAATCAGAAATGCATATTGGCTGTGAGCAAGAAGGTTCTTGGCTCGCCCGGGTGATAACGGTAGCCGCTCTTGTTGATTGAGGCGACGTAATCGGTATCAAACAGGTTGTAAACGTTTAGCTGGAAGTCGAGATTGCGATTAACTCGATACCCCAGTTTGGCATCGGCGACCCAGTAACCTTCGGTAAACGCTGGCGTTCCCACCGCGCCGTCTGAACCTTTATGCATACTGCCGATATAGCGTGCGCCCGCGCCAACAGAGATATCGTCGGTTGCCTGATATTGGCTCCATAAGGTGAAGGCGTGCTCCGGGGTATACGGCAGCGATGAGGAACCATCCTGGGCAACATCTTTGCCGTTTTTGATGGTTGCTTTTTGCTGGGTATAGCCGCCAATCACCTGCCACGCGGGAGTGATATTCCCGGCCACGGATATCTCATAGCCTTCGACGCGTTTCTTACCGTATTGCGAGTAAGTTCCGTCATCATTTTGCTCAACTTCATTTTCGATATCAGTGCGGAACAGCGCGGCGGTGAGCAACAGACGTTTATCCAGAACCTGCCATTTGGTGCCAATCTCGCTGGTGTTGGCTTTTTGCGGTTTAAAATCGGTGCGGTTGGCACTGTTACCGCTGCCAGACTGCGCAAGGGCGAAGTTGTTGCCGCCCGGAGGCTGCTGGGAAACGGCATAGTTAATATAGACATTGCCGTTTTCCGTCAGGTGATACAGCGCCCCGGCTTTCCAGTTCATCAGATTGCCCGACTTGGCGGTGTCGACGGTGGTGACCGGAGAACCTTTTGCCACACCAGTTGGGCAGGTGATGGCACCGCGTCCGCTGCCGCCGCAGGCGGTGGCACTGTCATATTCAGTATGATAATTATCCAGACGGATCCCGCCGTTCAGCTCAAAATCACGGGTGATTTGCAGCGTATCAAAGGCGTAAATTGCGAAGGTA	*
NC_000913.3_sliding:3988801-3989800	914	666	1000	-	papilloma.Late_protein_L2.Papillomavirus_panthera7600:BK066462	414	495	525	32.2	9.13e-04	8M1I9M2I19M2I5M4D37M	QLPQLADQLAALGESDLLFALSQHAVAFAQSQLHQQDRKWPRLPDYFAIGRTTALALHTVSGQKILYPQDREISEVLLQLPEL	TCGGCTTTTTGGCTCTCCTGGGCTTTTTGCAATGCCGTCAGTTGGTTAGCCAGGGCATCGCTGGTGGCGGTCTGATTGACGGCCTGTTGTTTACCCCAGCCATACAAACCGATGCCCGCCGCCAGAGCAATAGCGATAGCCACCGCGCTGAGAATCAATGCGGTATTGTTCTTACTCTTTTTTTCTGTTGCGACAGGTTGTGACGTGGTGTCCACGGCCTCCCTGGTCTCTTCAACCACGGCGGAGGTTTTTTCTTGTTCCGTCATTATGGCTTCCTGTTATGAGAGTTATTGTAATGCCCGTAAAAGCGCATCGTTGTCAGCGTTATCGGCGACCTTAATGTCTTGCCAGCCCAGTTCCCGGGCGAGTTTCGCCAAACGCTCACTGACGACCAATAGTCGACAGTGTAGTAACCAGTGCTCACGATACCATTGTGGGATCAGCGACCAGAGTTGCTGCAACATTTCACCGCTGGTAACAACGACCATCGTCACCTCGCGGGCTTGCCAGCGCATCGCTTCTTCTGCACCATCGTAATGGATTGCGCATCGTTGATAACATTCACAAAAAGTGACCTCAGCACCGCGCGCCGTCAGGGTATCCCCAATTAGCTCACGACCACCATTGCCACGTAATATCAGCGCACGTTTGCCCGCAATATTTTGTAATTCAGGTAATTGTAGCAAGACTTCGCTGATTTCCCGATCCTGCGGGTAGAGAATCTTCTGTCCACTTACGGTATGTAGTGCCAGTGCGGTGGTGCGTCCAATGGCGAAATAATCAGGTAGTCGGGGCCATTTACGATCTTGCTGATGCAGCTGTGATTGGGCAAAAGCAACCGCGTGTTGCGAGAGGGCAAACAACAGATCGCTCTCCCCCAGCGCTGCCAGTTGATCAGCAAGTTGCGGTAATTGTTGACCCGGAGAAAACTCAATCAGCGGAAAATGCCAGGCCACCTGCCCCAGTGTGCGCAGACGGCTCACTAACTCTTCTCCAGCGG	*
NC_000913.3_sliding:4339201-4340200	319	5	1000	-	papilloma.Late_protein_L1.Bandicoot_papillomatosis_carcinomatosis_virus_type_2:EU277647	237	339	506	28.0	5.45e-05	9M2I23M1D20M2I20M1D29M	DLLEKTSDRLHFDEAWYGYARFNPIYADHYAMRGEPGDHNGPTVFATHSTHKLLNALSQASYIHVREGRGAINFSRFNQAYMMHATTSPLYAICASNDVAVSMMD	TGCCGTCCATCATCGACACCGCCACGTCGTTGGATGCGCAGATGGCATACAGCGGGGAGGTGGTGGCATGCATCATGTAGGCCTGGTTGAAGCGGGAGAAGTTAATCGCCCCACGACCTTCACGTACATGAATATAAGAAGCCTGTGACAGCGCATTCAGCAGTTTGTGGGTGGAGTGGGTGGCGAAAACGGTAGGACCGTTGTGATCGCCAGGTTCGCCGCGCATGGCATAGTGATCGGCATAGATCGGGTTGAAACGTGCATAGCCGTACCAGGCTTCGTCAAAGTGCAGACGATCGGAGGTTTTTTCCAGCAGATCCTGCGCTTCTTTAGCGTTATAACACACGCCGTCATAGGTGCAGTTGGTCACCACGCAGTAAGACGGTTTTTGCCCGGCTTTGTCTTTGGTCAGCGGGCTTTCACTGATTTTCTTCTGCAAGGTTTCAGGTTGCATTTCCTGCGGATAGATTGGCCCGATAATGCCGTAGCGGTTGCGGCTTGGCACCATATAGACCGGTTTCGCGCCTGTCAGCATCAAACCTTGTTCGATGGATTTATGGCAGTTACGGTCAACGACCACGACATCGTTATCGGTCATGCAAGCCTGCATGATGGTGCGGTTAGAGCCGGAAGTACCGACGACTACCGACCAGGAGCGATCGGCACCAAATACGCGTGCGGCATATTTTTCGCTTTCGCCAAATGCGCCAGTATGGTCAAGCAAAGAACCGAGGGAAGTTCGTTCGATGCCCATGTCGGTGCGGAACAGATTTTCACCATAGTAGTCATGGTAGAAACGTCCGGCGGGTGTTTTGGTAAAACCAACGCCGCCCTGGTGGCCTGGCGCTGCCCAGGAATATTCATGGATGTCACTATATTTCATCAGCGCGCTGAACAGTGGCGGCAACAGCTGCTGGCGGTAGCGGGTCATCGCGGCAACGGCGCGTCCGGCGATAAAGTCGGCGGTATCTTCCAGAATCCAGGCGAATTCATCGACAAG	*
```
