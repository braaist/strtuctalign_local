#include "pdb.h"
#include <math.h>
#include <assert.h>
#include <float.h>

struct pdb_info {
	char Name;
	unsigned int all_m;
	unsigned int n;
	struct atom *all_atoms_dna;
	struct atom *atoms_prot;
	struct atom *atoms_wat;
	unsigned int n_pairs;
	char **pairs;
	int **compl_pairs;
	unsigned int *compl_list;
};

struct best_score {
	unsigned int *list_P1, *list_C11, *list_OP11, *list_OP21, *list_P2, *list_C12, *list_OP12, *list_OP22;
	unsigned int dna1_chain1_n, dna1_chain2_n, dna2_chain1_n, dna2_chain2_n;
	unsigned int i_max_measure, j_max_measure;
	unsigned int compl1, compl2;
	unsigned int n_first_chain, m_first_chain;
	unsigned int pair1, pair2;
	double S_max;
	struct atom *atoms_dna1, *atoms_dna2;
	struct atom *dna1_chain1, *dna1_chain2, *dna2_chain1, *dna2_chain2;
	unsigned int m1, m2;
};

void deleteScore(struct best_score* score) {
	free(score->list_P1);
	free(score->list_C11);
	free(score->list_OP11);
	free(score->list_OP21);
	free(score->list_P2);
	free(score->list_C12);
	free(score->list_OP12);
	free(score->list_OP22);
	deleteAtomlist(score->atoms_dna1, score->m1);
	deleteAtomlist(score->atoms_dna2, score->m2);
	deleteAtomlist(score->dna1_chain1, score->dna1_chain1_n);
	deleteAtomlist(score->dna1_chain2, score->dna1_chain2_n);
	deleteAtomlist(score->dna2_chain1, score->dna2_chain1_n);
	deleteAtomlist(score->dna2_chain2, score->dna2_chain2_n);
}

void findBestScore(const struct pdb_info* pdb1, const struct pdb_info* pdb2, struct best_score* result) {
	result->S_max = -DBL_MAX;
	/*** MAIN FOR CYCLE ***/
	struct atom *atoms_prot_CA1 = NULL;
	struct atom *atoms_prot_C1 = NULL;
	struct atom *atoms_prot_i1;
	struct atom *atoms_prot_j1, *atoms_prot_j2;
	struct atom *C_atoms_prot_i1 = NULL, *C_atoms_prot_j2 = NULL;
	for (unsigned int pair1=1; pair1<=pdb1->n_pairs; pair1++)
	{
		struct atom *atoms_dna1 = NULL, *dna1_chain1 = NULL, *dna1_chain2 = NULL;
		unsigned int m1, dna1_chain1_n=0, dna1_chain2_n=0;
		SelectChain(pdb1->all_atoms_dna, pdb1->all_m, &dna1_chain1, &dna1_chain1_n, pdb1->pairs[pair1][1]);
		SelectChain(pdb1->all_atoms_dna, pdb1->all_m, &dna1_chain2, &dna1_chain2_n, pdb1->pairs[pair1][2]);
		atomlistmerge(&atoms_dna1, &m1, dna1_chain1, dna1_chain1_n, dna1_chain2, dna1_chain2_n);


		/* Make lists of P, C1', OP1 atoms
		   of dna and CA atoms of protein */
		unsigned int *list_P1, *list_C11, *list_OP11, *list_OP21, *list_CA1, *list_C1;
		unsigned int n_P1, n_C11, n_OP11, n_OP21, n_CA1, n_C1;
		getAtomsNumbers(atoms_dna1, m1, &list_P1, &n_P1, "P");
		// pdb.c function, get only indexes of atoms
		getAtomsNumbers(atoms_dna1, m1, &list_C11, &n_C11, "C1'");
		getAtomsNumbers(atoms_dna1, m1, &list_OP11, &n_OP11, "OP1");
		getAtomsNumbers(atoms_dna1, m1, &list_OP21, &n_OP21, "OP2");
		getAtomsNumbers(pdb1->atoms_prot, pdb1->n, &list_CA1, &n_CA1, "CA");
		getAtomsNumbers(pdb1->atoms_prot, pdb1->n, &list_C1, &n_C1, "C");
		//assert(n_C1 == n_CA1);
		correctC1_P(atoms_dna1, &list_C11, &n_C11, list_P1, &n_P1);
		correctC_CA(pdb1->atoms_prot, &list_C1, &n_C1, &list_CA1, &n_CA1);
		//for (i=1; i<=n_C11; i++) printf("%s.%c\n", atoms_dna1[list_C11[i]].ResNumber, atoms_dna1[list_C11[i]].Chain);
		//corrects list_C1 to use rule list_P[i] and list_C1[i+1] are in the same nucleotide

		atoms_prot_CA1 = (struct atom *)malloc( sizeof(struct atom)*(n_CA1+1) );
		for (unsigned int i=1; i<=n_CA1; i++) {
			atomcpy(&atoms_prot_CA1[i], pdb1->atoms_prot[list_CA1[i]]);
			//pdb.c function. Copy all properties of atom
		}
		atoms_prot_C1 = (struct atom *)malloc( sizeof(struct atom)*(n_C1+1) );
		for (unsigned int i=1; i<=n_C1; i++) {
			atomcpy(&atoms_prot_C1[i], pdb1->atoms_prot[list_C1[i]]);
		}
		// printf("Done atoms:\tP %u\tC1 %u\tOP1 %u\tOP2 %u\tCA %u\tC %u\n",n_P1, n_C11, n_OP11, n_OP21, n_CA1, n_C1);
		for (unsigned int pair2=1; pair2<=pdb2->n_pairs; pair2++)
		{
			//printf("\nCYCLE %c:%c vs %c:%c\n", pdb1->pairs[pair1][1], pdb1->pairs[pair1][2], pdb2->pairs[pair2][1], pdb2->pairs[pair2][2]);
			struct atom *atoms_dna2 = NULL, *dna2_chain1, *dna2_chain2;
			unsigned int m2, dna2_chain1_n=0, dna2_chain2_n=0;
			SelectChain(pdb2->all_atoms_dna, pdb2->all_m, &dna2_chain1, &dna2_chain1_n, pdb2->pairs[pair2][1]);
			SelectChain(pdb2->all_atoms_dna, pdb2->all_m, &dna2_chain2, &dna2_chain2_n, pdb2->pairs[pair2][2]);
			atomlistmerge(&atoms_dna2, &m2, dna2_chain1, dna2_chain1_n, dna2_chain2, dna2_chain2_n);



			unsigned int *list_P2, *list_C12, *list_OP12, *list_OP22, *list_CA2, *list_C2;
			unsigned int n_P2, n_C12, n_OP12, n_OP22, n_CA2, n_C2;

			getAtomsNumbers(atoms_dna2, m2, &list_P2, &n_P2, "P");
			getAtomsNumbers(atoms_dna2, m2, &list_C12, &n_C12, "C1'");
			getAtomsNumbers(atoms_dna2, m2, &list_OP12, &n_OP12, "OP1");
			getAtomsNumbers(atoms_dna2, m2, &list_OP22, &n_OP22, "OP2");
			getAtomsNumbers(pdb2->atoms_prot, pdb2->n, &list_CA2, &n_CA2, "CA");
			getAtomsNumbers(pdb2->atoms_prot, pdb2->n, &list_C2, &n_C2, "C");
			//assert(n_C2 == n_CA2);
			correctC1_P(atoms_dna2, &list_C12, &n_C12, list_P2, &n_P2);
			correctC_CA(pdb2->atoms_prot, &list_C2, &n_C2, &list_CA2, &n_CA2);

			struct atom *atoms_prot_CA2 = (struct atom *)malloc( sizeof(struct atom)*(n_CA2+1) );
			for (unsigned int i=1; i<=n_CA2; i++) {
				atomcpy(&atoms_prot_CA2[i], pdb2->atoms_prot[list_CA2[i]]);
			}
			struct atom *atoms_prot_C2 = (struct atom *)malloc( sizeof(struct atom)*(n_C2+1) );
			for (unsigned int i=1; i<=n_C2; i++) {
				atomcpy(&atoms_prot_C2[i], pdb2->atoms_prot[list_C2[i]]);
			}
			// printf("Done atoms:\tP %u\tC1 %u\tOP1 %u\tOP2 %u\tCA %u\tC %u\n",n_P2, n_C12, n_OP12, n_OP22, n_CA2, n_C2);

			/* Done making lists */

			/* Create array of measures for
			   all pairs of dna P atoms (list_measure).
			   dim(list_measure) = n_P2 x n_P1 */
			double **list_measure;
			unsigned int n_hit;
			unsigned int **list_hit;

			list_measure = (double **)malloc( (n_P2+1)*sizeof(double *)*(n_P1+1) );
			for (unsigned int i=1; i<=n_P1; i++)  list_measure[i] = (double *)malloc( sizeof(double)*(n_P2+1) );
			//for (i=1; i<=n_P1; i++) for (j=1; j<=n_P2;j++)  list_measure[i][j] = (unsigned int *)malloc( sizeof(unsigned int) );
			// Enable if measure function returns unsigned int
			for (unsigned int i=1; i<=n_P1; i++){
				for (unsigned int j=1; j<=n_P2; j++){
					list_measure[i][j] = 0;
					if (i > n_OP11 || i > n_OP21 || i+1 > n_C11)continue;
					if (j > n_OP12 || j > n_OP22 || j+1 > n_C12)continue;
					ChangeSystem(atoms_prot_CA1, n_CA1, &atoms_prot_i1, atoms_dna1[list_P1[i]], atoms_dna1[list_C11[i+1]], atoms_dna1[list_OP11[i]], atoms_dna1[list_OP21[i]], 'E');
					// pdb.c function. Change the coordinate system of protein with given nucleotide
					ChangeSystem(atoms_prot_C1, n_C1, &C_atoms_prot_i1, atoms_dna1[list_P1[i]], atoms_dna1[list_C11[i+1]], atoms_dna1[list_OP11[i]], atoms_dna1[list_OP21[i]], 'E');
					ChangeSystem(atoms_prot_CA2, n_CA2, &atoms_prot_j2, atoms_dna2[list_P2[j]], atoms_dna2[list_C12[j+1]], atoms_dna2[list_OP12[j]], atoms_dna2[list_OP22[j]], 'F');
					ChangeSystem(atoms_prot_C2, n_C2, &C_atoms_prot_j2, atoms_dna2[list_P2[j]], atoms_dna2[list_C12[j+1]], atoms_dna2[list_OP12[j]], atoms_dna2[list_OP22[j]], 'F');
					BidirectionalHit(atoms_prot_i1, C_atoms_prot_i1, n_CA1, atoms_prot_j2, C_atoms_prot_j2, n_CA2, &list_hit, &n_hit);
					// pdb.c function.
					Measure2_p(&(list_measure[i][j]), list_hit, n_hit, atoms_prot_i1, atoms_prot_j2);
					// pdb.c function.
					//printf("Measure: %f  i: %u j: %u n_P1: %u  n_P2: %u \n",(list_measure[i][j]), i,j, n_P1, n_P2);
					// Enable in test mode
					deleteAtomlist(atoms_prot_i1, n_CA1);
					deleteAtomlist(C_atoms_prot_i1, n_C1);
					deleteAtomlist(atoms_prot_j2, n_CA2);
					deleteAtomlist(C_atoms_prot_j2, n_C2);
					for(int i=0;i<=n_CA1;i++)
						free(list_hit[i]);
					free(list_hit);
				}
			}



			/* Done creation of array of measures */


			/* Start working with diagonals */

			unsigned int i_max, j_max, i_start, j_start, i_max_measure, j_max_measure;
			unsigned int n_first_chain, m_first_chain;
			unsigned int compl1=0, compl2=0;
			double S_max;
			struct atom *atoms_dna_P1 = NULL;
			struct atom *atoms_dna_P2 = NULL;
			struct atom *atoms_dna_C11 = NULL;
			struct atom *atoms_dna_C12 = NULL;



			// print measure-table
			/*
			printf("  ");
			  for (int j=1; j<=n_P2; j++) printf("%4d", j);
			  puts("");
			  for (int i=n_P1; i>=1; i--){
			  printf("%2d", i);
			  for (int j=1; j<=n_P2; j++){
			  printf("%4.0f", list_measure[i][j]>0 ? list_measure[i][j] : 0);
			  }
			  puts("");
			  }
			  printf("  ");
			  for (int j=1; j<=n_P2; j++) printf("%4d", j);
			  puts("");
			  */

			i_max_measure=0;
			j_max_measure=0;
			atoms_dna_P1 = (struct atom *)malloc( sizeof(struct atom)*(n_P1+1) );
			for (unsigned int i=1; i<=n_P1; i++) {
				atomcpy(&atoms_dna_P1[i], atoms_dna1[list_P1[i]]);
			}
			atoms_dna_P2 = (struct atom *)malloc( sizeof(struct atom)*(n_P2+1) );
			for (unsigned int i=1; i<=n_P2; i++) {
				atomcpy(&atoms_dna_P2[i], atoms_dna2[list_P2[i]]);
			}
			atoms_dna_C12 = (struct atom *)malloc( sizeof(struct atom)*(n_C12+1) );
			for (unsigned int i=1; i<=n_C12; i++) {
				atomcpy(&atoms_dna_C12[i], atoms_dna2[list_C12[i]]);
			}
			atoms_dna_C11 = (struct atom *)malloc( sizeof(struct atom)*(n_C11+1) );
			for (unsigned int i=1; i<=n_C11; i++) {
				atomcpy(&atoms_dna_C11[i], atoms_dna1[list_C11[i]]);
			}
			//find_compl(atoms_dna1, list_P1, list_C11, list_OP11, list_OP21, atoms_dna_P1, n_P1, &compl1, &n_first_chain);
			(run_find_compl(atoms_dna_C11, n_C11, &compl1, &n_first_chain, pdb1->compl_pairs[pair1]));
			//find_compl(atoms_dna2, list_P2, list_C12, list_OP12, list_OP22, atoms_dna_P2, n_P2, &compl2, &m_first_chain);
			(run_find_compl(atoms_dna_C12, n_C12, &compl2, &m_first_chain, pdb2->compl_pairs[pair2]));
			BestDiag(list_measure, n_P1, n_P2, &S_max, &i_max, &j_max, &i_start, &j_start, &i_max_measure, &j_max_measure,
					atoms_dna1, list_P1, atoms_dna2, list_P2, compl1, compl2, n_first_chain, m_first_chain);
			// pdb.c function
			/* Done diagonal search */

			if (S_max > result->S_max)
			{
				result->S_max = S_max;
				atomlistcpy(&result->dna1_chain1, dna1_chain1, dna1_chain1_n);
				result->dna1_chain1_n = dna1_chain1_n;
				atomlistcpy(&result->dna1_chain2, dna1_chain2, dna1_chain2_n);
				result->dna1_chain2_n = dna1_chain2_n;
				atomlistcpy(&result->dna2_chain1, dna2_chain1, dna2_chain1_n);
				result->dna2_chain1_n = dna2_chain1_n;
				atomlistcpy(&result->dna2_chain2, dna2_chain2, dna2_chain2_n);
				result->dna2_chain2_n = dna2_chain2_n;

				atomlistcpy(&result->atoms_dna1, atoms_dna1, m1);
				result->m1 = m1;
				result->list_P1 = (unsigned int *)malloc(sizeof(unsigned int)*(n_P1+1));
				result->list_C11 = (unsigned int *)malloc(sizeof(unsigned int)*(n_C11+1));
				result->list_OP11 = (unsigned int *)malloc(sizeof(unsigned int)*(n_OP11+1));
				result->list_OP21 = (unsigned int *)malloc(sizeof(unsigned int)*(n_OP21+1));
				memcpy(result->list_P1, list_P1, sizeof(unsigned int)*(n_P1+1));
				memcpy(result->list_C11, list_C11, sizeof(unsigned int)*(n_C11+1));
				memcpy(result->list_OP11, list_OP11, sizeof(unsigned int)*(n_OP11+1));
				memcpy(result->list_OP21, list_OP21, sizeof(unsigned int)*(n_OP21+1));

				atomlistcpy(&result->atoms_dna2, atoms_dna2, m2);
				result->m2 = m2;
				result->list_P2 = (unsigned int *)malloc(sizeof(unsigned int)*(n_P2+1));
				result->list_C12 = (unsigned int *)malloc(sizeof(unsigned int)*(n_C12+1));
				result->list_OP12 = (unsigned int *)malloc(sizeof(unsigned int)*(n_OP12+1));
				result->list_OP22 = (unsigned int *)malloc(sizeof(unsigned int)*(n_OP22+1));
				memcpy(result->list_P2, list_P2, sizeof(unsigned int)*(n_P2+1));
				memcpy(result->list_C12, list_C12, sizeof(unsigned int)*(n_C12+1));
				memcpy(result->list_OP12, list_OP12, sizeof(unsigned int)*(n_OP12+1));
				memcpy(result->list_OP22, list_OP22, sizeof(unsigned int)*(n_OP22+1));
				result->i_max_measure = i_max_measure;
				result->j_max_measure = j_max_measure;
				result->compl1 = compl1;
				result->compl2 = compl2;
				result->n_first_chain = n_first_chain;
				result->m_first_chain = m_first_chain;
				result->pair1 = pair1;
				result->pair2 = pair2;
			}
			free(atoms_prot_CA2);
			free(atoms_prot_C2);
			free(atoms_dna_P1);
			free(atoms_dna_P2);
			free(atoms_dna_C11);
			free(atoms_dna_C12);
			free(dna2_chain1);
			free(dna2_chain2);
			free(atoms_dna2);
			for (int i=1;i<=n_P1;i++) {
				free(list_measure[i]);
			}
			free(list_measure);
		}
		free(atoms_prot_CA1);
		free(atoms_prot_C1);
			free(dna1_chain1);
			free(dna1_chain2);
			free(atoms_dna1);
	}
	/*** END OF MAIN CYCLE ***/
	assert(result->S_max > -DBL_MAX);
}


void reverseAtoms(struct atom* atoms, size_t n) {
	for (size_t i=0; i*2<n; i++) {
		struct atom tmp = atoms[i];
		atoms[i] = atoms[n-1-i];
		atoms[n-1-i] = tmp;
	}
}


int main  (int argc, char **argv)
{
	puts("       ______               __  ___   ___             ");
	puts("      / __/ /_______ ______/ /_/ _ | / (_)__ ____     ");
	puts("     _\\ \\/ __/ __/ // / __/ __/ __ |/ / / _ `/ _ \\    ");
	puts("    /___/\\__/_/  \\_,_/\\__/\\__/_/ |_/_/_/\\_, /_//_/    ");
	puts("                                       /___/          ");
	puts("\n");
	puts("*********************************************************");
	puts("* Program StructAlign inputs a pair of complexes of DNA *");
	puts("* double helix with proteins and outputs an alignment   *");
	puts("* of DNA chains corresponding to the best spatial fit   *");
	puts("* of the protein chains. StructAlign is developed and   *");
	puts("* maintained by a team from A.N. Belozersky Institute   *");
	puts("* of Physico-Chemical Biology MSU and Faculty of        *");
	puts("* Bioengineering and Bioinformatics MSU.                *");
	puts("*********************************************************");
	puts("\n");
	/* Checking command line,
	   throw exception if not complete */
	printf("%s %s %s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

	if (argc != 9)
	{
		printf("\nUsage: %s -i <input file> -o <output file> -f <file name for score> -g <gk file>", argv[0]);
		printf("\nExample: %s -i input.txt -o out.txt -f results.txt -g out.gk\n\n", argv[0]);
		exit (1);
	}

  FILE * fp;
  char *outfile, *inpfile, *gkfile, *max_score_filename;
  char** infiles;
  max_score_filename = (char *)malloc( sizeof(char)*(strlen(argv[6])+1) );
  inpfile = (char *)malloc( sizeof(char)*(strlen(argv[2])+1) );
  outfile = (char *)malloc( sizeof(char)*(strlen(argv[4])+1) );
  gkfile = (char *)malloc( sizeof(char)*(strlen(argv[8])+1) );
  infiles = malloc(sizeof(char*) * 12);
  char* chain_names = malloc(12);
	char** starts = malloc(sizeof(char*) * 12);
	char** ends = malloc(sizeof(char*) * 12);

  sscanf(argv[2], "%s", inpfile);
  sscanf(argv[4], "%s", outfile);
  sscanf(argv[6], "%s", max_score_filename);
  sscanf(argv[8], "%s", gkfile);

  fp = fopen(inpfile, "r");
  if (fp == NULL)
    exit(1);

  char str1[10], str2[10], str3[10], str4[10];
  unsigned long num_pdbs = 0;
  while (fscanf(fp, "%s %s %s %s\n", str1, str2, str3, str4) != EOF){
    infiles[num_pdbs] = malloc(strlen(str1)+1);
    strcpy(infiles[num_pdbs], str1);
    sscanf(str2, "%c", &chain_names[num_pdbs]);
		starts[num_pdbs] = malloc(strlen(str3)+1);
		strcpy(starts[num_pdbs], str3);
		ends[num_pdbs] = malloc(strlen(str4)+1);
		strcpy(ends[num_pdbs], str4);
    num_pdbs++;
		}
	/* Done reading arguments */

	/* DECLARATION AND ASSIGNMENT
	   OF SOME VARIABLES */
	/*** For server ***/
	FILE *max_score;
	max_score = fopen(max_score_filename, "w");
	if (max_score == NULL) {
		perror(max_score_filename);
		exit(1);
	}

	struct pdb_info* pdbs = malloc(sizeof(struct pdb_info) * num_pdbs);
	for (size_t pdb_i = 0; pdb_i < num_pdbs; pdb_i++) {
		unsigned int maxnumber = 256; //primary number of atoms in each groups (DNA, protein, water)
		FILE *flow_out;
		/* m - number of DNA atoms
		 * n - number of protein atoms
		 * w - number of water atoms */

		/* Reading PDB file */
		struct atomname *list;
		char *prot_chains = (char *)malloc( sizeof(char)*(2+1) );
		char *dna_chains = (char *)malloc( sizeof(char)*(2+1) );
		unsigned int w=0;

		printf("Reading %zu PDBs file...\n\n", pdb_i+1);
		struct pdb_info *cur_pdb = &pdbs[pdb_i];
		readerPDB(infiles[pdb_i], &cur_pdb->all_m, &cur_pdb->n, &w, maxnumber, &cur_pdb->all_atoms_dna, &dna_chains, &cur_pdb->atoms_prot, &prot_chains, &cur_pdb->atoms_wat, &list);
		// pdb.c function, read PDB and put protein, DNA and water atoms in three different arrays
		printf("\n...done; %d atoms of dna, %d atoms of protein, %d atoms of water\n", cur_pdb->all_m, cur_pdb->n, w);
		if (cur_pdb->all_m == 0)
		{	fprintf(max_score, "Error\n%zu structure has no DNA!", pdb_i);
			exit(1); }
		if (chain_names[pdb_i]=='@')
			chain_names[pdb_i] = prot_chains[1];
		else
			if ( inArray(chain_names[pdb_i], prot_chains)==0 )
			{ fprintf(max_score, "Error\nChain %c is not a protein chain in first structure", chain_names[pdb_i]);
				exit(1); }
		SelectChain(cur_pdb->atoms_prot, cur_pdb->n, &cur_pdb->atoms_prot, &cur_pdb->n, chain_names[pdb_i]);
		SelectRange(cur_pdb->atoms_prot, cur_pdb->n, &cur_pdb->atoms_prot, &cur_pdb->n, starts[pdb_i], ends[pdb_i]);
		// pdb.c function
		printf("Atoms in selected chain: %u\n", cur_pdb->n);
		// print to stdout number of protein atoms in selected chain

		dna_chains[0] = ' ';
		/*** For server ***/

		printf("DNA chains identified in %s : %s;\n", str1, dna_chains);
		/* Done reading */

		/*** 3DNA block ***/
		printf("%s\n", "Running find_pair from x3dna suite:");
		run_3dna(infiles[pdb_i], &cur_pdb->compl_list, &cur_pdb->compl_pairs, &cur_pdb->pairs, &cur_pdb->n_pairs, max_score_filename, '_', '_');
		assert(cur_pdb->n_pairs > 0);
		printf("%s\n", "find_pair completed succesfully!\n");
	}
	/*** 3DNA block end ***/
	printf("SEARCHING FOR CENTROID\n");
	size_t centroid = -1;
	double max_min_score = -1;
	for (size_t pdb1_i = 0; pdb1_i < num_pdbs; pdb1_i++) {
		printf("TRYING PDB %zu/%zu\n", pdb1_i+1, num_pdbs);
		struct pdb_info* pdb1 = &pdbs[pdb1_i];
		double min_score = 1e70;
		for (size_t pdb2_i = 0; pdb2_i < num_pdbs; pdb2_i++) {
			if (pdb1_i == pdb2_i) continue;
			struct pdb_info* pdb2 = &pdbs[pdb2_i];
			struct best_score probe;
			findBestScore(pdb1, pdb2, &probe);
			if (probe.S_max < min_score) {
				min_score = probe.S_max;
			}
			deleteScore(&probe);
		}
		if (min_score > max_min_score) {
			max_min_score = min_score;
			centroid = pdb1_i;
		}
	}
	printf("USE PDB %zu/%zu AS CENTROID (MIN SCORE IS %f)\n", centroid+1, num_pdbs, max_min_score);


	char max_chain = 'A';
	struct atom** chains = malloc(sizeof(struct atom*) * num_pdbs * 3);
	unsigned int * n_chains = malloc(sizeof(unsigned int) * num_pdbs * 3);
	for (size_t i = 0; i < num_pdbs; i++) {
		//if (i == centroid) continue;
		struct best_score probe;
		findBestScore(&pdbs[i], &pdbs[centroid], &probe);
		unsigned int dna_n11, dna_n12, dna_n21, dna_n22; //sequence length
		char *dna_seq11, *dna_seq12, *dna_seq21, *dna_seq22;
		char *dna_num_seq11, *dna_num_seq12, *dna_num_seq21, *dna_num_seq22;
		int dna1_chain1_start, dna1_chain2_start, dna2_chain1_start, dna2_chain2_start;

		Seq(probe.dna1_chain1, probe.dna1_chain1_n, &dna_seq11, &dna_num_seq11, &dna_n11, max_score);
		sscanf(probe.dna1_chain1[1].ResNumber, "%d", &dna1_chain1_start);

		Seq(probe.dna1_chain2, probe.dna1_chain2_n, &dna_seq12, &dna_num_seq12, &dna_n12, max_score);
		sscanf(probe.dna1_chain2[1].ResNumber, "%d", &dna1_chain2_start);

		Seq(probe.dna2_chain1, probe.dna2_chain1_n, &dna_seq21, &dna_num_seq21, &dna_n21, max_score);
		sscanf(probe.dna2_chain1[1].ResNumber, "%d", &dna2_chain1_start);

		Seq(probe.dna2_chain2, probe.dna2_chain2_n, &dna_seq22, &dna_num_seq22, &dna_n22, max_score);
		sscanf(probe.dna2_chain2[1].ResNumber, "%d", &dna2_chain2_start);

		if (probe.S_max == 0)
			fprintf(max_score, "Warning! The score is zero! It seems you have specified wrong parameters.\n");

		unsigned int i_max_measure_compl, j_max_measure_compl;
		unsigned int is_reverse1, is_reverse2;

		//printf("n_first_chain=%u i_max_measure=%u best_compl1=%u\n", best_n_first_chain, best_i_max_measure,  best_compl1);
		is_reverse1 = ((probe.i_max_measure > probe.n_first_chain) ? 1 : 0);
		is_reverse2 = ((probe.j_max_measure > probe.m_first_chain) ? 1 : 0);
		i_max_measure_compl = probe.compl1-probe.i_max_measure+1;
		j_max_measure_compl = probe.compl2-probe.j_max_measure+1;
		//printf("i_max_measure=%u i_max_measure_compl=%u\n", best_i_max_measure, i_max_measure_compl);
		//printf("i_max_measure_compl=%s\n", best_atoms_dna1[best_list_P1[i_max_measure_compl]].ResNumber);
		if (i_max_measure_compl > dna_n11+dna_n12)
		{
			j_max_measure_compl = j_max_measure_compl-(i_max_measure_compl-dna_n11-dna_n12);
			i_max_measure_compl = dna_n11+dna_n12;
		}
		if (j_max_measure_compl > dna_n21+dna_n22)
		{
			i_max_measure_compl = i_max_measure_compl-(j_max_measure_compl-dna_n21-dna_n22);
			j_max_measure_compl = dna_n21+dna_n22;
		}

		int i_max_measure_compl_num = 0;
		int j_max_measure_compl_num = 0;
		char *i_max_measure_compl_str, *j_max_measure_compl_str;
		i_max_measure_compl_str = (char *)malloc( sizeof(char)*7 );
		j_max_measure_compl_str = (char *)malloc( sizeof(char)*7 );
		sscanf(probe.atoms_dna1[probe.list_P1[i_max_measure_compl]].ResNumber, "%d", &i_max_measure_compl_num);
		sscanf(probe.atoms_dna2[probe.list_P2[j_max_measure_compl]].ResNumber, "%d", &j_max_measure_compl_num);

		if (i_max_measure_compl==probe.i_max_measure)
		{
			sscanf(probe.atoms_dna1[probe.list_P1[i_max_measure_compl+1]].ResNumber, "%d", &i_max_measure_compl_num);
			i_max_measure_compl_num = i_max_measure_compl_num - 1;
			//puts("Here1");
		}
		if (j_max_measure_compl==probe.j_max_measure)
		{
			sscanf(probe.atoms_dna2[probe.list_P2[j_max_measure_compl+1]].ResNumber, "%d", &j_max_measure_compl_num);
			j_max_measure_compl_num = j_max_measure_compl_num - 1;
			//puts("Here2");
		}


		int start;
		start = (is_reverse1 == 1) ? dna1_chain1_start : dna1_chain2_start;
		if (i_max_measure_compl_num < start)
		{
			//printf("ires_compl=%d start=%d\n", i_max_measure_compl_num, start);
			j_max_measure_compl_num = j_max_measure_compl_num + start - i_max_measure_compl_num;
			i_max_measure_compl_num = start;
			//printf("ires_compl=%d start=%d\n", i_max_measure_compl_num, start;
		}
		start = (is_reverse2 == 1) ? dna2_chain1_start : dna2_chain2_start;
		if (j_max_measure_compl_num < start)
		{
			//printf("jres_compl=%d start=%d\n", j_max_measure_compl_num, start);
			i_max_measure_compl = i_max_measure_compl_num + start - j_max_measure_compl_num;
			j_max_measure_compl = start;
			//printf("jres_compl=%d start=%d\n", j_max_measure_compl_num, start);
		}

		sprintf(i_max_measure_compl_str, "%d", i_max_measure_compl_num);
		sprintf(j_max_measure_compl_str, "%d", j_max_measure_compl_num);

		fprintf(max_score, "%c %c %c %c %c %c %s %s %s %s %u %u\n%lg\n%s %s\n%s %s\n%s %s\n%s %s\n", chain_names[i], chain_names[centroid], pdbs[i].pairs[probe.pair1][1], pdbs[i].pairs[probe.pair1][2], pdbs[centroid].pairs[probe.pair2][1], pdbs[centroid].pairs[probe.pair2][2], probe.atoms_dna1[probe.list_P1[probe.i_max_measure]].ResNumber, probe.atoms_dna2[probe.list_P2[probe.j_max_measure]].ResNumber, i_max_measure_compl_str, j_max_measure_compl_str, is_reverse1, is_reverse2, probe.S_max, dna_seq11, dna_num_seq11, dna_seq12, dna_num_seq12, dna_seq21, dna_num_seq21, dna_seq22, dna_num_seq22);

		struct atom *atoms_dna_i1 = NULL;
		struct atom *atoms_dna_i2 = NULL;
		struct atom *atoms_dna_i1_tmp = NULL;
		struct atom *atoms_dna_i2_tmp = NULL;
		/* Change the system to i and j coordinates, write the alignment to file */

		struct atom *atoms_prot_i1;
		struct atom *atoms_prot_i1_tmp;
		struct atom *atoms_prot_j1, *atoms_prot_j2;

		struct coordsystem invertion = ChangeSystem(pdbs[centroid].atoms_prot, pdbs[centroid].n, &atoms_prot_j2, probe.atoms_dna2[probe.list_P2[probe.j_max_measure]], probe.atoms_dna2[probe.list_C12[probe.j_max_measure+1]], probe.atoms_dna2[probe.list_OP12[probe.j_max_measure]], probe.atoms_dna2[probe.list_OP22[probe.j_max_measure]], 'F');
		struct coordsystem prot_system = ChangeSystem(pdbs[i].atoms_prot, pdbs[i].n, &atoms_prot_i1_tmp, probe.atoms_dna1[probe.list_P1[probe.i_max_measure]], probe.atoms_dna1[probe.list_C11[probe.i_max_measure+1]], probe.atoms_dna1[probe.list_OP11[probe.i_max_measure]], probe.atoms_dna1[probe.list_OP21[probe.i_max_measure]], 'E');
		ChangeSystemR(atoms_prot_i1_tmp, pdbs[i].n, &atoms_prot_i1, probe.atoms_dna2[probe.list_P2[probe.j_max_measure]], max_chain++, invertion);
		struct coordsystem dna1_system = ChangeSystem(probe.dna1_chain1, probe.dna1_chain1_n, &atoms_dna_i1_tmp, probe.atoms_dna1[probe.list_P1[probe.i_max_measure]], probe.atoms_dna1[probe.list_C11[probe.i_max_measure+1]], probe.atoms_dna1[probe.list_OP11[probe.i_max_measure]], probe.atoms_dna1[probe.list_OP21[probe.i_max_measure]], 'A');
		ChangeSystemR(atoms_dna_i1_tmp, probe.dna1_chain1_n, &atoms_dna_i1, probe.atoms_dna2[probe.list_P2[probe.j_max_measure]], max_chain++, invertion);
		struct coordsystem dna2_system = ChangeSystem(probe.dna1_chain2, probe.dna1_chain2_n, &atoms_dna_i2_tmp, probe.atoms_dna1[probe.list_P1[probe.i_max_measure]], probe.atoms_dna1[probe.list_C11[probe.i_max_measure+1]], probe.atoms_dna1[probe.list_OP11[probe.i_max_measure]], probe.atoms_dna1[probe.list_OP21[probe.i_max_measure]], 'B');
		ChangeSystemR(atoms_dna_i2_tmp, probe.dna1_chain2_n, &atoms_dna_i2, probe.atoms_dna2[probe.list_P2[probe.j_max_measure]], max_chain++, invertion);
		// pdb.c function
		chains[i*3+0] = malloc(sizeof(struct atom) * (pdbs[i].n + 1));
		memcpy(chains[i*3+0], atoms_prot_i1, sizeof(struct atom) * (pdbs[i].n + 1));
		n_chains[i*3+0] = pdbs[i].n;
		chains[i*3+1] = malloc(sizeof(struct atom) * (probe.dna1_chain1_n + 1));
		memcpy(chains[i*3+1], atoms_dna_i1, sizeof(struct atom) * (probe.dna1_chain1_n+1));
		n_chains[i*3+1] = probe.dna1_chain1_n;
		chains[i*3+2] = malloc(sizeof(struct atom) * (probe.dna1_chain2_n + 1));
		memcpy(chains[i*3+2], atoms_dna_i2, sizeof(struct atom) * (probe.dna1_chain2_n + 1));
		n_chains[i*3+2] = probe.dna1_chain2_n;
		/*if (invertion.ort_X.X*prot_system.ort_X.X < 0) {
			reverseAtoms(chains[i*3+0]+1, pdbs[i].n);
		}
		if (invertion.ort_X.X*dna1_system.ort_X.X < 0) {
			reverseAtoms(chains[i*3+1]+1, probe.dna1_chain1_n);
		}
		if (invertion.ort_X.X*dna2_system.ort_X.X < 0) {
			reverseAtoms(chains[i*3+2]+1, probe.dna1_chain2_n);
		}*/
	}
	fclose(max_score);
	createPDB(outfile, outfile);
	for(int i=0;i<num_pdbs*3;i++) {
		writetoPDB(outfile, chains[i], n_chains[i]);
	}
	endPDB(outfile);
	writegk(gkfile, outfile, infiles, num_pdbs, chains, n_chains, num_pdbs*3, centroid);
	return 0;
}
/* End main */
