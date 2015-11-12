This is the old version of our GeneQuery tool -both stand-alone used for bootstrapping and other applications, and web-server. 

Also includes C code used to make two necessary binary files: fisher_test_right and normal_pval_left

Compiling commands: 

gcc -std=c99 -O2 normal.c normal_pval_left.c -o normal_pval_left -lm

gcc -std=c99 -O2 fisher.c   fisher_test_right.c  -o fisher_test_right  -lm

Perl requirements: CGI, File::Temp, HTML::Table

Images and database itself are of substantial size and thus are not deposited here (but are available upon request). 
