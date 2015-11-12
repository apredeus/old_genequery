This is the old version of our GeneQuery tool -both stand-alone used for bootstrapping and other applications, and web-server. 

1. The repository includes C code used to make two necessary binary files: fisher_test_right and normal_pval_left

Compiling commands: 

gcc -std=c99 -O2 normal.c normal_pval_left.c -o normal_pval_left -lm

gcc -std=c99 -O2 fisher.c   fisher_test_right.c  -o fisher_test_right  -lm

2. Perl requirements: CGI, File::Temp, HTML::Table

3. Necessary data. The datais  of substantial size and thus are not deposited here (but are available upon request). 

Includes:
	- images for 
	- databases in GSEXXX_GPLYYY<tab>Module#<tab>Gene_entrez_id format (1 line per gene). 
	- 3-column annotation for each GPL that's on the list - tab-separated probe ID, gene symbol, and Entrez ID. 

