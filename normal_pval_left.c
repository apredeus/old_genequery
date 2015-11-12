#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>


// Wrapper module for one-sided p-value calculation for normal distribution -- A.P.
// Done using R functions used to calculate p-value (written in C) 
// Added functions to evaluate average reference logp form Fisher's test (mu) and standard deviation (sigma). 

double main_pval(double fisher_logp, uint32_t db, uint32_t module_size); 
double get_mu(uint32_t db, uint32_t module_size);
double get_sigma(uint32_t db, uint32_t module_size); 
double pnorm5(double x, double mu, double sigma, uint32_t lower_tail, uint32_t log_p);
void pnorm_both(double x, double *cum, double *ccum, uint32_t i_tail, uint32_t log_p);
#define MAXLINELEN 131072

int main(int argc, char** argv) {
  FILE* test_file;
  char buf[MAXLINELEN];
  char* bufptr;
  
  double   fisher_logp; 
  uint32_t db;
  uint32_t module_size;
  char s1[MAXLINELEN]; 
  char s2[MAXLINELEN]; 
  char s3[MAXLINELEN]; 
  char s4[MAXLINELEN]; 
  char s5[MAXLINELEN]; 

  if (argc != 4) {
  main_std_help:
    printf(
"One-sided p-value for normally distributed variables.\n"
"Usage: normal_pval_left [filename] [db] [module size]\n\n"
"[db] is an integer identifier:\n" 
"	mm_2K database = \"1\"\n"
"	mm_4K database = \"2\"\n"
"	hs_2K database = \"3\"\n"
"	hs_4K database = \"4\"\n" 
"each line of the filename is expected to have log p-value, calculated by Fisher's exact test,\n"  
"in the 6th position of the tab-delimited file. This can be changed if necessary. \n\n"
"Written by Alexander Predeus using functions by Ross Ihaka, R Core team, and the R foundation (see pnorm.c for more) .\n"
	   );
    return 1;
  }

// Read the file and get values 

  test_file = fopen(argv[1], "r");
  db = atoi(argv[2]);
  module_size = atoi(argv[3]);

  if (!test_file) {
    printf("Error: Unable to open file.\n");
    return 2;
  }

  buf[MAXLINELEN - 1] = ' ';
  while (fgets(buf, MAXLINELEN, test_file)) {
    if (!buf[MAXLINELEN - 1]) {
      printf("Error: Excessively long line in input file.\n");
      fclose(test_file);
      return 3;
    }
    bufptr = buf;
    while ((*bufptr == ' ') || (*bufptr == '\t')) {
      bufptr++;
    }
    if (*bufptr < ' ') {
      continue;
    }
    if (sscanf(bufptr, "%s %s %s %s %s %lf", s1, s2, s3, s4, s5, &fisher_logp) == 6 ){
      printf("%s\t%s\t%s\t%s\t%s\t%.3f\t%e\t%.3f\n",s1,s2,s3,s4,s5,fisher_logp,main_pval(fisher_logp, db, module_size),log10(main_pval(fisher_logp, db, module_size)));
    } else {
      // skip improperly formatted line
      continue;
    }
  }
  if (!feof(test_file)) {
    printf("Error: File read failure.\n");
    fclose(test_file);
    return 4;
  }
  fclose(test_file);
  return 0;
}
