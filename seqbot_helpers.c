#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "seqbot_helpers.h"

/* Return the melting temperature of sequence, or -1 if the sequence is invalid.
 * The melting temperature formula is given in the handout.
 * An invalid sequence is a sequence of length 0, or a sequence that contains
 * characters other than 'A', 'C', 'G', 'T'.
 */
int calculate_melting_temperature(char *sequence, int sequence_length)
{
	if (sequence_length < 0){
		return -1;
	
	}
	int n_A = 0;
	int n_T = 0;
	int n_C = 0;
	int n_G = 0;

	int i;
	for (i = 0; i < sequence_length; i++){
		if(sequence[i] == 'A'){
			n_A += 1;
		} 
		else if (sequence[i] == 'C'){
			n_C += 1;
		
		}
		else if (sequence[i] == 'G'){
			n_G += 1;

		}
		else if (sequence[i] == 'T'){
			n_T += 1;
		
		}
		else{
			return -1;
		}
	}
	
	int temp = (n_A + n_T) * 2 + (n_C + n_G) * 4;
        return temp;
}

/* Prints the instructions to make a molecule from sequence.
 * If an invalid character is found in sequence print
 * "INVALID SEQUENCE" and return immediately
 */
void print_instructions(char *sequence, int sequence_length)
{
	if(sequence_length < 0){
		printf("INVALID SEQUENCE");
		exit(1);
	}
	printf("START\n");
	char prev = sequence[0];
	int prev_n = 0;
	int i;
	for(i = 0; i < sequence_length; i++){
		if (sequence[i] == 'A'){
			if (prev != 'A'){
				printf("WRITE %c %d\n", prev, prev_n);
				prev = 'A';
				prev_n = 0;
			}
			if (i == sequence_length - 1){
				printf("WRITE %c %d\n", prev, prev_n + 1);
			
			}

			else{
				prev_n += 1;
			}
	  	}	
		else if (sequence[i] == 'C'){
			if (prev != 'C'){
				printf("WRITE %c %d\n", prev, prev_n);
				prev = 'C';
				prev_n = 0;
			}
			if (i == sequence_length - 1){
				printf("WRITE %c %d\n", prev, prev_n + 1);
			}
			else{
				prev_n += 1;
			}
		
		}
		else if (sequence[i] == 'G'){
			if(prev != 'G'){
				printf("WRITE %c %d\n", prev, prev_n);
				prev = 'G';
				prev_n = 0;
			}
			if(i == sequence_length - 1){
				printf("WRITE %c %d\n", prev, prev_n + 1);
			}
			else{
				prev_n += 1;
			}
		}
	   	else if(sequence[i] == 'T'){
	   		if(prev != 'T'){
				printf("WRITE %c %d\n", prev, prev_n);
				prev = 'T';
				// 1 since we have already seen 1 T nucleotide
				prev_n = 0;
			}
			if (i == sequence_length - 1){
				printf("WRITE %c %d\n", prev, prev_n + 1);
			}
			else{
				prev_n += 1;
			}
		}
	   
	   
	   else{
	   	printf("INVALID SEQUENCE");
		exit(1); 
	   }
		
      }
	int temp = calculate_melting_temperature(sequence, sequence_length);
	printf("SET_TEMPERATURE %d\n", temp);
        printf("END\n");
}


/* Print to standard output all of the sequences of length k.
 * The format of the output is "<length> <sequence> 0" to 
 * correspond to the input format required by generate_molecules_from_file()
 * 
 * Reminder: you are not allowed to use string functions in these files.
 */

void comb_helper(int index, int k, char *molec);

void generate_all_molecules(int k)
{

	if(k > 0){
		char molec[k];	
		comb_helper(0, k, molec);
	}
	else{
		printf("INVALID SEQUENCE");
		exit(1);
	}
}

/*helper for getting all combination of nucleotides*/

void comb_helper(int index, int k, char *molec){

	char nuc[4] = {'A', 'T', 'C', 'G'};
	if ((index + 1) == k){
		int i;
		for (i = 0; i < 4; i ++){
			molec[index] = nuc[i];
			//printf("%d %s 0\n", k, molec);
			printf("%d ", k);
			int s;
			for (s = 0; s < k; s++){
				printf("%c", molec[s]);
			}	
			printf(" 0\n");
		}		
	}

	else{
		int j;
		for(j = 0; j < 4; j++){
			molec[index] = nuc[j];
			comb_helper(index + 1, k, molec);
		}
	}




}

/*helper to complement the strand*/
void helper_complement(char *string, int nuc_size){
	char complem_str[nuc_size];
	for (int j = 0; j < nuc_size; j++){
		if (string[j] == 'A'){
			complem_str[j] = 'T';}
		
		else if (string[j] == 'T'){
			complem_str[j] = 'A';}
	
	
		else if (string[j] == 'C'){

			complem_str[j] = 'G';}
								

		else if (string[j] == 'G'){
			complem_str[j] = 'C';}

	}


	print_instructions(complem_str, nuc_size);

	

}

/*Helper to reverse the strand*/

void helper_reverse(char *string, int nuc_size){

	char rev_str[nuc_size];

	for(int j = 0; j < nuc_size; j++){
		rev_str[j] = string[nuc_size - j - 1]; 
	}

	print_instructions(rev_str, nuc_size);

}

/*Helper to reverse and compelement the strand*/
void helper_rev_complem(char *string, int nuc_size){

	char rev_str[nuc_size];
	for(int j = 0; j < nuc_size; j++){
	
		rev_str[j] = string[nuc_size - j - 1];
	}
	helper_complement(rev_str, nuc_size);
	
	

}

	
	/* Print the instructions for each of the sequences found in filename according
 * to the mode provided.
 * filename contains one sequence per line, and the format of each line is
 * "<length> <sequence> <mode>" where 
 *     - <length> is the number of characters in the sequence 
 *     - <sequence> is the array of DNA characters
 *     - <mode> is either 0, 1, 2, or 3 indicating how the <sequence> should 
 *              be modified before printing the instrutions. The modes have the 
 *              following meanings:
 *         - 0  - print instructions for sequence (unmodified)
 *         - 1  - print instructions for the the complement of sequence
 *         - 2  - print instructions for the reverse of sequence
 *         - 3  - print instructions for sequence where it is complemented 
 *                and reversed.
 * 
 * Error checking: If any of the following error conditions occur, the function
 * immediately prints "INVALID SEQUENCE" to standard output and exits the 
 * program with a exit code of 1.
 *  - length does not match the number of characters in sequence
 *  - length is not a positive number
 *  - sequence contains at least one invalid character
 *  - mode is not a number between 0 and 3 inclusive
 * 
 * You do not need to verify that length or mode are actually integer numbers,
 * only that they are in the correct range. It is recommended that you use a 
 * fscanf to read the numbers and fgetc to read the sequence characters.
 */
void generate_molecules_from_file(char* filename)
{

	int nuc_size;
	int mode;
	FILE *fp = fopen(filename, "r");
	
	if (fp == NULL){
		fprintf(stderr, "Error Opening File\n");
		exit(1);
	}

	while(fscanf(fp, "%d ", &nuc_size) == 1){

		if (nuc_size <= 0){
		
			printf("INVALID SEQUENCE");
			exit(1);
		}
		char string[nuc_size];
		int i;
		for (i = 0; i < nuc_size; i++){
			
			string[i] = fgetc(fp);

			// captures if size of nuc is not right, or if the characters are not right
			if ((string[i] != 'A') && (string[i] != 'C') && (string[i] != 'T') && (string[i] != 'G')){
			
				printf("INVALID SEQUENCE");
				exit(1);
			
			}

		}

		// sequence longer than intended
		char extra_c = fgetc(fp);
		if ( (extra_c == 'A') | (extra_c == 'C') | (extra_c == 'T') | (extra_c == 'G')){
			
			printf("INVALID SEQUENCE");
			exit(1);
		}	



		fscanf(fp, "%d ", &mode);

		if(mode == 0){
			print_instructions(string, nuc_size);
		
		}
		else if(mode == 1){

			helper_complement(string, nuc_size);
			
			}
		else if (mode == 2){
			helper_reverse(string, nuc_size);
		}
		else if (mode == 3){
			helper_rev_complem(string, nuc_size);
		
		}
		else{
			printf("INVALID SEQUENCE");
			exit(1);
		
		}
		
	
	}


	fclose(fp);
    }

    




