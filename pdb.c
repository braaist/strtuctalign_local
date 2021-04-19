#include "pdb.h"
#include "nts.c"
#include <math.h>
#include <assert.h>
#define TRUE 1
#define FALSE 0
#define EPSILON 0.000001


#include <time.h>
#include <stdlib.h>


/********************************
*		ADDITIONAL WORK WITH READ PDB
*********************************/

unsigned int SelectChain(struct atom * atoms_from, unsigned int n_from, struct atom ** atoms_to, unsigned int * n_to, char chain)
  {
  unsigned int i, j=0;

  (*atoms_to) = (struct atom *)malloc( sizeof(struct atom)*(n_from+1) );
  /*printf("OK1\n");*/ // Enable in test mode
  for (i=1; i<=n_from; i++) {
    if (atoms_from[i].Chain == chain) {
	  j++;
	  (*atoms_to)[j] = atoms_from[i];
	  }
	}
  /*printf("OK\n");*/ // Enable in test mode
  *n_to = j;

  return 0;
}

unsigned int SelectRange(struct atom * atoms_from, unsigned int n_from, struct atom ** atoms_to, unsigned int * n_to, char *start, char *end)
  {
  unsigned int i, j=0;
  double buff, s, e;

  (*atoms_to) = (struct atom *)malloc( sizeof(struct atom)*(n_from+1) );
  /*printf("OK1\n");*/ // Enable in test mode

  if (strcmp(start, "zero") == 0)
	s = -INFINITY;
  else
	sscanf(start, "%lg", &s);
  if (strcmp(end, "inf") == 0 )
	e = INFINITY;
  else
	sscanf(end, "%lg", &e);

  for (i=1; i<=n_from; i++) {
    sscanf(atoms_from[i].ResNumber, "%lg", &buff);
    //printf("%lg %s %s\t", buff, start, end);
  if (buff >= s && buff <= e) {
    j++;
    (*atoms_to)[j] = atoms_from[i];
  }
/*
    if (strcmp(start, "zero") == 0)
    {
    	if (strcmp(end, "inf") != 0)
    	{
    		sscanf(end, "%d", &e);
    		if (buff <= e) {
			  j++;
			  (*atoms_to)[j] = atoms_from[i];
		}
	}
	else
	{
		j++;
		(*atoms_to)[j] = atoms_from[i];
	}
    }
    else
    {
    	sscanf(start, "%d", &s);
    	//printf("!!!%d %d!!!", buff, s);
	if (strcmp(end, "inf") == 0)
	{
    		if (buff >= s) {
			  j++;
			  (*atoms_to)[j] = atoms_from[i];
		}
	}
	else
	{
		sscanf(end, "%d", &e);
		if (buff >= s && buff <= e) {
		  j++;
		  (*atoms_to)[j] = atoms_from[i];
		}
	}
    }*/
  }

  /*printf("OK\n");*/ // Enable in test mode
  *n_to = j;

  return 0;
}

unsigned int getAtomsNumbers(struct atom * atoms_from, unsigned int n, unsigned int ** list_to, unsigned int * n_to, char * atoms_type){
  unsigned int i, k;

  *list_to = (unsigned int *)malloc( sizeof(unsigned int)*(n+1) );
  k=1;
  for (i=1; i<=n; i++){
    if (strcmp(atoms_from[i].AtomName, atoms_type) == 0){(*list_to)[k] = i; k++;}
  }
  *n_to = k-1;
  return 0;
 }

unsigned int correctC_CA(struct atom *atoms_from, unsigned int **list_C, unsigned int *n_C, unsigned int **list_CA, unsigned int *n_CA) {
	for (unsigned int i=1; i <= *n_C; i++) {
		if (i > *n_CA) {
			*n_C = *n_CA;
		}
    		if ( strcmp(atoms_from[(*list_C)[i]].ResNumber, atoms_from[(*list_CA)[i]].ResNumber) != 0 ) {
    			if (i + 1 <= *n_CA && strcmp(atoms_from[(*list_C)[i]].ResNumber, atoms_from[(*list_CA)[i+1]].ResNumber) == 0 ) {
				for (int j = i; j < *n_CA; j++) {
					(*list_CA)[j] = (*list_CA)[j+1];
				}
				(*n_CA)--;
			}
			if (i - 1 > 0 && strcmp(atoms_from[(*list_C)[i]].ResNumber, atoms_from[(*list_CA)[i-1]].ResNumber) == 0 ) {
				for (int j = i-1; j < *n_C; j++) {
					(*list_C)[j] = (*list_C)[j+1];
				}
				(*n_C)--;
				i--;
			}
		}
	}
	*n_CA = *n_C;
  return 0;
}
unsigned int correctC1_P(struct atom *atoms_from, unsigned int **list_C, unsigned int *n_C, unsigned int *list_P, unsigned int *n_P)
  {
  unsigned int i, j;

  for (j=1; j<=(*n_P); j++)
  {
    //printf("P: %s; C1: %s\n", atoms_from[list_P[j]].ResNumber, atoms_from[(*list_C)[j+1]].ResNumber);
    if ( j == (*n_C) ) //if number of P-atoms is greater than C1'-atoms
    {
    	(*n_P) -= 1;
    	return 0;
    }

    if ( strcmp(atoms_from[list_P[j]].ResNumber, atoms_from[(*list_C)[j+1]].ResNumber) != 0 )
    {

      if ( strcmp(atoms_from[list_P[j]].ResNumber, atoms_from[(*list_C)[j+2]].ResNumber) == 0 )
      {
        (*n_C) -= 1;
        for (i=j+1; i<=(*n_C); i++)
          {(*list_C)[i] = (*list_C)[i+1];}
        //printf("P: %s; C1: %s\n", atoms_from[list_P[j]].ResNumber, atoms_from[(*list_C)[j+1]].ResNumber); //prints if smth changed
      }

      if ( strcmp(atoms_from[list_P[j]].ResNumber, atoms_from[(*list_C)[j]].ResNumber) == 0 )
      {
        (*n_C) += 1;
        (*list_C) = (unsigned int *)realloc( (*list_C), sizeof(unsigned int)*((*n_C)+1) );
        for (i=(*n_C); i>j; i--)
          {(*list_C)[i] = (*list_C)[i-1];}
        //printf("P: %s; C1: %s\n", atoms_from[list_P[j]].ResNumber, atoms_from[(*list_C)[j+1]].ResNumber); //prints if smth changed
      }

    }
  }
  return 0;
}


/********************************
*   For server
********************************/
int inArray(char c, char *string)
{
	unsigned int i;
	unsigned int leng = 1;

	while (string[leng] != '\0')
		leng++;

	for (i=1; i<leng; i++)
	{
		if (string[i]==c)
			return TRUE;
	}
	return FALSE;
}

void Seq(struct atom *atoms, unsigned int n, char **seq, char **num_seq, unsigned int *m,  FILE *max_score)
{
	unsigned int i;
	unsigned int *list_C;
	char *res;
	char nt;
	res = (char *)malloc( sizeof(char)*4 );
	getAtomsNumbers(atoms, n, &list_C, &(*m), "C1'");
	(*seq) = (char *)malloc( sizeof(char)*((*m)+1) );
	(*num_seq) = (char *)malloc( sizeof(char)*(*m)*5);
	**num_seq = 0;
	for (i=1; i<=(*m); i++)
	{
		strcpy(res, atoms[list_C[i]].ResType);
		//fprintf(max_score, "%c:%s.%s.%s\n", atoms[list_C[i]].Chain, atoms[list_C[i]].AtomName, atoms[list_C[i]].ResType, atoms[list_C[i]].ResNumber);
		if ( strcmp(res, "DT")==0 || strcmp(res, "DA")==0 || strcmp(res, "DG")==0 || strcmp(res, "DC")==0 )
		{
			(*seq)[i-1] = res[1];
		}
		else
		{
			if ( strcmp(res, "U")==0 || strcmp(res, "A")==0 || strcmp(res, "G")==0 || strcmp(res, "C")==0 )
			{
				(*seq)[i-1] = res[0];
			}
			else
			{
				if ( convertNT( res, &nt) == 0 )
				{
					(*seq)[i-1] = nt;
					fprintf(max_score, "Warning! There is modified residue %s`%s in chain %c of %s. Reading it as %c.\n", res, atoms[list_C[i]].ResNumber, atoms[list_C[i]].Chain, atoms[list_C[i]].PDBCode, nt);
				}
				else
				{
					fprintf(max_score, "Warning! The program doesn't know residue %s\n! Please, report us about it.\n", res);
					(*seq)[i-1] = '?';
				}
			}
		}
		strcat( (*num_seq), atoms[list_C[i]].ResNumber );
		if ( i<(*m) )
			strcat( (*num_seq), "," );
	}
	(*seq)[(*m)] = '\0';
}

/********************************
*		MEASURES FOR HIT COUNT
*********************************/

double Measure1(unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2)
  {
  unsigned int i,j;
  double measure = 0.0;
  struct geomvector vector;
  for (i=1; i<=n_hit; i++) {
    for (j=i; j<=n_hit; j++){
	  vector = fromto( location(atoms_prot_1[list_hit[i][0]]), location(atoms_prot_2[list_hit[j][1]]) );
	  if (scalar(vector, vector) > 0.0) measure += 1/sqrt( scalar(vector, vector) );
	  /*printf("m1:\t%f\tsqrt:\t%f\n", measure, sqrt( scalar(vector, vector) ));*/ // Enable in test mode
	}
  }
  return measure;
}

double Measure2(unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2)
  {
  unsigned int i,j;
  double measure = 0.0;
  struct geomvector vector;
  for (i=1; i<=n_hit; i++) {
	  vector = fromto( location(atoms_prot_1[list_hit[i][0]]), location(atoms_prot_2[list_hit[i][1]]) );
	  measure += 1.5 - sqrt( scalar(vector, vector) );
	  printf("measure:\t%f\tsqrt:\t%f\n", measure, sqrt( scalar(vector, vector) ));
  }
  printf("I return measure: %f\n", measure);
  return measure;
}

double Measure3(unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2)
  {
	  /*printf("m3:\t%u\n", n_hit);*/ // Enable in test mode
  return n_hit;
}

unsigned int Measure1_p(double *measure, unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2)
  {
  unsigned int i,j;
  double m = 0.0;
  struct geomvector vector;
  for (i=1; i<=n_hit; i++) {
    for (j=i; j<=n_hit; j++){
	  vector = fromto( location(atoms_prot_1[list_hit[i][0]]), location(atoms_prot_2[list_hit[j][1]]) );
	  if (scalar(vector, vector) > 0.0) m += 1/sqrt( scalar(vector, vector) );
	  /*printf("m1:\t%f\tsqrt:\t%f\n", measure, sqrt( scalar(vector, vector) ));*/ // Enable in test mode
	}
  }
  * measure = m;
  return 0;
}

double Measure2_p(double *measure, unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2)
  {
  unsigned int i,j;
  double m = 0.0;
  struct geomvector vector;
  for (i=1; i<=n_hit; i++) {
	  vector = fromto( location(atoms_prot_1[list_hit[i][0]]), location(atoms_prot_2[list_hit[i][1]]) );
	  m += ( 4.5 - sqrt(scalar(vector, vector))>0 ? 4.5 - sqrt(scalar(vector, vector)) : 0 );
	  /*printf("measure:\t%f\tsqrt:\t%f\n", m, sqrt( scalar(vector, vector) ));*/ // Enable in test mode
  }
  *measure = m-31.0193;
  /*printf("I return measure: %f\n", *measure);*/ // Enable in test mode
  return 0;
}

double Measure3_p(double *measure, unsigned int ** list_hit, unsigned int n_hit, struct atom * atoms_prot_1, struct atom * atoms_prot_2)
  {
	  /*printf("m3:\t%u\n", n_hit);*/ // Enable in test mode
  *measure = n_hit;
  return 0;
}

/********************************
*                            3DNA
********************************/
unsigned int run_3dna(char *pdb_name, unsigned int **compl, int ***compl_pairs, char ***pairs, unsigned int *len, unsigned int server, char *max_score_filename, char priority1, char priority2)
{
	FILE *fp, *out_file;
	extern FILE *popen();
	char *command;
	char c[102];
	char chain1, chain2, flag;
	unsigned int pairs_max=2;
	unsigned int i, j, n, a, b, count;
	int res1, res2;
	command = (char *)malloc(sizeof(char)*(strlen(pdb_name)+150));
	(*compl) = (unsigned int *)malloc(sizeof(unsigned int)*(pairs_max+1));
	(*compl_pairs) = (int **)malloc(sizeof(int *)*(pairs_max+1));
	(*pairs) = (char **)malloc(sizeof(char *)*(pairs_max+1));
	for (j=0; j<=pairs_max; j++)
	{
		(*pairs)[j] = (char *)malloc(sizeof(char)*3);
		(*compl_pairs)[j] = (int *)malloc(sizeof(int)*3);
	}
	if (server == 0)
	{
		sprintf(command, "find_pair %s out 2>/dev/null", pdb_name);
		fp = popen(command, "r");
		if (fp == NULL)
		{
			printf("Failed to run find_pair on %s\n", pdb_name);
			exit(1);
		}
		pclose(fp);
		out_file = fopen("out", "r");
		if (out_file == NULL)
		{
			perror("find_pair failed");
			exit(1);
		}
	}
	else
	{
    sprintf(command, "find_pair %s find_pair.output", pdb_name);
		fp = popen(command, "r");
		if (fp == NULL)
		{
			FILE *max_score;
			max_score = popen(max_score_filename, "w");
				fprintf(max_score, "Error\nFailed to run find_pair on %s\n", pdb_name);
				exit(1);
		}
		pclose(fp);
		char outname[146];
		sprintf(outname, "find_pair.output");
		out_file = fopen(outname, "r");
		if (out_file == NULL)
		{
			FILE *max_score;
			max_score = popen(max_score_filename, "w");
			fprintf(max_score, "Error\nfind_pair failed");
			exit(1);
		}
	}
	fgets (c, 102, out_file);
	fgets (c, 102, out_file);
	fgets (c, 102, out_file);
	fgets (c, 102, out_file);
	sscanf(c, "%u", &n);
	if (n==0)
	{
		perror("There is no double stranded DNA!");
		exit(1);
	}
	fgets (c, 102, out_file);

	flag = 'x';
	count = 0;
	chain1 = '@';
	chain2 = '@';
	(*pairs)[count][1] = chain1;
	(*pairs)[count][2] = chain2;
	int foundPriority = 0;
	char tmp;
	for (i=1; i<=n; i++)
	{
		fgets (c, 102, out_file);
		sscanf(c, "%5u%5u%*u #%*u %c %*4c>%c%*[:.]%d%*19c%*[:.]%d%*c:%c", &a, &b, &tmp, &chain1, &res1, &res2,  &chain2);
		if (chain1 != (*pairs)[count][1] && chain2 != (*pairs)[count][2]){
			flag = 'x';
		}
		if (!foundPriority && (chain1 == priority1 && chain2 == priority2 || chain2 == priority1 && chain1 == priority2)) {
			flag = 'x';
			foundPriority = 1;
		}

		if (flag == 'x')
		{
			count++;
			if (count > pairs_max)
			{
				pairs_max = pairs_max * 2;
				(*compl) = (unsigned int *)realloc((*compl), (pairs_max+1)*sizeof(unsigned int));
				(*compl_pairs) = (int **)realloc((*compl_pairs), sizeof(int *)*(pairs_max+1));
				(*pairs) = (char **)realloc((*pairs), sizeof(char *)*(pairs_max+1));
				for (j=pairs_max; j>pairs_max/2; j--)
				{
					(*pairs)[j] = (char *)malloc(sizeof(char)*3);
					(*compl_pairs)[j] = (int *)malloc(sizeof(int)*3);
				}
			}
			(*compl)[count] = b-a; // unnessesary
			(*compl_pairs)[count][1] = res1;
			(*compl_pairs)[count][2] = res2;
			(*pairs)[count][1] = chain1;
			(*pairs)[count][2] = chain2;
			(*len) = count;
			n --;
			flag = ' ';
		}
	}
	return 0;
}


/********************************
*		BIDIRECTIONAL HIT
*********************************/


unsigned int BidirectionalHit( struct atom * atoms_i, struct atom *C_atoms_i, unsigned int n_i, struct atom * atoms_j, struct atom *C_atoms_j, unsigned int n_j, unsigned int *** list, unsigned int * n_hit){
  unsigned int i,j,k, curr, m, *list_j, *list_i, ja, ma;
  double min, dist;
  struct geomvector vector, vector1, vector2;

  /*free (*list);*/
  //printf("H1\n"); // Enable in test mode
  *list = (unsigned int **)malloc( 2*sizeof(unsigned int *)*(n_i+2) );
  for (i=0; i<n_i+1; i++) (*list)[i] = (unsigned int *)malloc( sizeof(unsigned int)*2 );
  //printf("H2\n"); // Enable in test mode

  list_i = (unsigned int *)malloc(sizeof(unsigned int)*(n_i+1));
  list_j = (unsigned int *)malloc(sizeof(unsigned int)*(n_j+1));
  //printf("H3\n"); // Enable in test mode

  m=0;
  for (i=1; i<=n_i; i++){
      min = 500.0;
      list_i[i] = 0;
      for (j=1; j<=n_j; j++){
	  vector = fromto( location(atoms_i[i]), location(atoms_j[j]) );
	  dist = sqrt( scalar(vector, vector) );
          if (dist - min < 0) {
          	vector1 = fromto( location(atoms_i[i]), location(C_atoms_i[i]) );
          	vector2 = fromto( location(atoms_j[j]), location(C_atoms_j[j]) );
          	if ( scalar(vector1, vector2) > 0 )
		{
			min = dist;
			list_i[i] = j;
		}
	  }
      }
      ma = 2;
      for (ja=1; ja<=m; ja++) {if (list_i[i]==list_j[ja]) ma=1;
      //printf("ma: %u\tlist_j[ja]:%u\tja: %u\n",ma, list_j[ja], ja); // Enable in test mode
      }
	  if (ma>1){
	    m++;
	    list_j[m] = list_i[i];
		//printf(" I add   m: %u\tlist_i[i]: %u\n",m, list_i[i]); // Enable in test mode
	    }
	//printf("\ti: %u list_i[i]: %u\tm: %u\tlist_j[m]: %u\n", i, list_i[i],m,list_j[m]); // Enable in test mode
  }
  //printf("H4\n"); // Enable in test mode

  k = 0;
  for (j=1; j<=m; j++){
      min = 500.0;
      for (i=1; i<=n_i; i++) {
	    vector = fromto( location(atoms_j[list_j[j]]), location(atoms_i[i]) );
	    dist = sqrt( scalar(vector, vector) );
        if (dist - min < 0) {
        	vector1 = fromto( location(atoms_i[i]), location(C_atoms_i[i]) );
          	vector2 = fromto( location(atoms_j[list_j[j]]), location(C_atoms_j[list_j[j]]) );
          	if ( scalar(vector2, vector1) > 0 )
		{
			min = dist;
			curr = i;
		}
	}
      }
    //printf("list_j[j]: %u\tlist_i[curr]: %u\tk: %u\t|\t",list_j[j],list_i[curr],k); // Enable in test mode
	if (list_j[j]==list_i[curr]){
	  //printf("I add: k: %u 0: %u 1: %u\n", k, curr, list_j[j]); // Enable in test mode
	  k++;
	  (*list)[k][0] = curr;
 	  (*list)[k][1] = list_j[j];
	  }
	/*if (min-90.0>0)printf("k: %u i: %u curr: %u min: %f \n", k, i, curr, min);*/ // Enable in test mode

	}
  //printf("H5\n");
  free(list_i);
  free(list_j);
  * n_hit = k;
  return 0;
  }

/********************************
*		CHANGE COORDINATE SYSTEM
*********************************/
// OP1 is an atom that has smaller scalar(P_OPx, P_C1)

struct coordsystem ChangeSystem(struct atom * atoms_from,
                          unsigned int n, struct atom ** atoms_to,
						  struct atom Op, struct atom Xp, struct atom Yp, struct atom Yp2,
						  char chainname)
{
  unsigned int i;
  struct geompoint Pp, Cpx, Cpy, Cpy2, Cpz, currp, zero;
  struct geompoint point3; //middle point between OP1 and OP2
  struct geomvector Cvx, Cvy, Cvy2, Cvz, currv;
  struct coordsystem dnares;
  zero.X = 0;
  zero.Y = 0;
  zero.Z = 0;

  //printf("L1\n");  // Enable in test mode
  Pp = location(Op);		// P coordinates
  Cpx = location(Xp);		// C1' coordinates
  Cvx = fromto(Pp, Cpx);
  if ( scalar(Cvx, Cvx) == 0 ) {
    perror("Cvx has zero length!!!");
    exit(1);
  }
  Cvx = ort(fromto(Pp, Cpx));	// P_C1' ort
  Cpy = location(Yp);		// OP1 in PDB
  Cpy2 = location(Yp2);		// OP2 in PDB

  point3.X = (Cpy.X + Cpy2.X)/2;
  point3.Y = (Cpy.Y + Cpy2.Y)/2;
  point3.Z = (Cpy.Z + Cpy2.Z)/2;
  Cvy = fromto(Pp, point3);
  if ( scalar(Cvy, Cvy) == 0 ){
    perror("Cvy has zero length!!!");
    exit(1);
  }
  Cvy = ort(Cvy);

  Cvz = vecproduct(Cvx, Cvy);
  if ( scalar(Cvz, Cvz) == 0 ) {
    perror("Cvz has zero length!!!");
    exit(1);
  }
  Cvz = ort(vecproduct(Cvx, Cvy));
  Cvy = vecproduct(Cvz, Cvx);
  if ( scalar(Cvy, Cvy) == 0 ) {
    perror("Cvy_vecproduct has zero length!!!");
    exit(1);
  }
  Cvy = ort(vecproduct(Cvz, Cvx));
  //printf("L2\n");  // Enable in test mode

  dnares.origin = zero;
  dnares.ort_X = Cvx;
  dnares.ort_Y = Cvy;
  dnares.ort_Z = Cvz;

  //printf("L3\n");  // Enable in test mode

  *atoms_to = (struct atom *)malloc( sizeof(struct atom)*(n+1) );
    //printf("L4\n");  // Enable in test mode

  for (i=1; i<=n; i++){
    atomcpy(&(*atoms_to)[i], atoms_from[i]);

	currp = location( (*atoms_to)[i] );
	currv = fromto( Pp, currp);
	currv = changesystem(currv, dnares, 1);

    (*atoms_to)[i].XCoord = currv.X;
    (*atoms_to)[i].YCoord = currv.Y;
    (*atoms_to)[i].ZCoord = currv.Z;
    (*atoms_to)[i].ChainRenamed = chainname;
  }
    //printf("L5\n");  // Enable in test mode

  return dnares;
  }

unsigned int ChangeSystemR(struct atom * atoms_from,
                          unsigned int n, struct atom ** atoms_to,
						  struct atom Op, char chainname, struct coordsystem dnares)
{
  unsigned int i;
  struct geompoint Pp, currp, zero;
  struct geomvector currv;
  zero.X = 0;
  zero.Y = 0;
  zero.Z = 0;

  Pp = location(Op);		// P coordinates


  *atoms_to = (struct atom *)malloc( sizeof(struct atom)*(n+1) );

  for (i=1; i<=n; i++){
    atomcpy(&(*atoms_to)[i], atoms_from[i]);

	currp = location( (*atoms_to)[i] );
	currv = changesystem(fromto(zero, currp), dnares, -1);
	currv = sumvec( currv, fromto(zero, Pp));

    (*atoms_to)[i].XCoord = currv.X;
    (*atoms_to)[i].YCoord = currv.Y;
    (*atoms_to)[i].ZCoord = currv.Z;
    (*atoms_to)[i].ChainRenamed = chainname;
  }

  return 0;
  }

/********************************
* 	BEST DIAGONAL FRAGMENT  *
*********************************/
void find_compl(struct atom *atoms1, unsigned int *list_P, unsigned int *list_C1, unsigned int *list_OP1, unsigned int *list_OP2, struct atom *atoms1_P, unsigned int n_P, unsigned int *compl, unsigned int *first_chain_length)
{
//atoms1 - all atoms
//atoms1_P - only P atoms

	unsigned i, k, c, count, k_min, n_P1, n_P2;
	double sum = 0;
	double min;
	double *sum_list;
	struct atom *atoms_dna_P = NULL; //for P-atoms in new coord system
	struct atom *atoms_first_chain = NULL; //for P-atoms in first chain
	struct atom *atoms_other_chain = NULL; //for P-atoms in not first chain

	sum_list = (double *)malloc( sizeof(double)*(n_P+1) );

	n_P1 = 0; n_P2 = 0;
	for (i=1; i<=n_P; i++)
	{
		//printf("%c.%s\n", atoms1_P[i].Chain, atoms1_P[i].ResNumber);
		if (atoms1_P[i].Chain == atoms1_P[1].Chain)
			{n_P1 += 1; }
		else {n_P2 += 1; }
	}

	(*first_chain_length) = n_P1;

	if (n_P2 == 0)
	{
		//printf(max_score, "Error\nThere is only 1 DNA chain!");
		perror("There is only 1 DNA chain!");
		exit(1);
	}

	atoms_first_chain = (struct atom *)malloc( sizeof(struct atom)*(n_P1+1) );
	atoms_other_chain = (struct atom *)malloc( sizeof(struct atom)*(n_P2+1) );
	sum_list = (double *)malloc( sizeof(double)*(n_P2+1) );

	c = 1;
	for (i=1; i<=n_P; i++)
	{
		if (atoms1_P[i].Chain == atoms1_P[1].Chain)
		{
			atomcpy(&atoms_first_chain[c], atoms1_P[i]);
			c += 1;
		}
		else
		{
			atomcpy(&atoms_other_chain[i-c+1], atoms1_P[i]);
		}
	}

	//for (k=n_P/2+2; k<=n_P+n_P/2; k++)
	for (k=2; k<=n_P; k++)
	{
		count = 0;
		sum = 0;
		for (i=1; i<=n_P1; i++)
		{

			ChangeSystem(atoms_other_chain, n_P2, &atoms_dna_P, atoms1[list_P[i]], atoms1[list_C1[i+1]], atoms1[list_OP1[i]], atoms1[list_OP2[i]], 'E');
			if ( (k-i>=1) && (k-i<=n_P2) )
			{
				//printf("basis1: P %s, C1 %s, OP1 %s, OP2 %s; relative: %s\n", atoms1[list_P[i]].ResNumber, atoms1[list_C1[i+1]].ResNumber, atoms1[list_OP1[i]].ResNumber, atoms1[list_OP2[i]].ResNumber, atoms_other_chain[k-i].ResNumber);
				count += 1;
				sum += pow(atoms_dna_P[k-i].XCoord-16.0526, 2)/2.1633 + pow(atoms_dna_P[k-i].YCoord+3.3171, 2)/3.1701 + pow(atoms_dna_P[k-i].ZCoord+8.5705, 2)/2.9010;

				//printf("basis2: P %s, C1 %s, OP1 %s, OP2 %s; relative: %s\n", atoms1[list_P[n_P1+k-i]].ResNumber, atoms1[list_C1[n_P1+k-i+1]].ResNumber, atoms1[list_OP1[n_P1+k-i]].ResNumber, atoms1[list_OP2[n_P1+k-i]].ResNumber, atoms_first_chain[i].ResNumber);
				ChangeSystem(atoms_first_chain, n_P1, &atoms_dna_P, atoms1[list_P[n_P1+k-i]], atoms1[list_C1[n_P1+k-i+1]], atoms1[list_OP1[n_P1+k-i]], atoms1[list_OP2[n_P1+k-i]], 'E');
				count += 1;
				sum += pow(atoms_dna_P[i].XCoord-16.0526, 2)/2.1633 + pow(atoms_dna_P[i].YCoord+3.3171, 2)/3.1701 + pow(atoms_dna_P[i].ZCoord+8.5705, 2)/2.9010;
			}
		sum_list[k] = sum/count;
		}
	}

	min = sum_list[2];
	k_min = 2;
	for (k=2; k<=n_P; k++)
	{
		if (sum_list[k] < min)
		{
			min = sum_list[k];
			k_min = k;
		}
	}

	k = k_min;
	//printf("k=%u\n", k);
	for (i=1; i<=n_P1; i++)
	{
		//printf("%u  %d\n", i, k-i);
		if ( (k-i>=1) && (k-i<=n_P2) && ( (strcmp(atoms1_P[i].ResType, "DA")==0) || (strcmp(atoms1_P[i].ResType, "DT")==0) || (strcmp(atoms1_P[i].ResType, "DG")==0) || (strcmp(atoms1_P[i].ResType, "DC")==0) ) )
		{
			if ( ( (strcmp(atoms1_P[i].ResType, "DA")==0 && strcmp(atoms_other_chain[k-i].ResType, "DA")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DG")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DC")!=0) ||
			(strcmp(atoms1_P[i].ResType, "DT")==0 && strcmp(atoms_other_chain[k-i].ResType, "DT")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DG")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DC")!=0) ||
			(strcmp(atoms1_P[i].ResType, "DG")==0 && strcmp(atoms_other_chain[k-i].ResType, "DA")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DG")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DT")!=0) ||
			(strcmp(atoms1_P[i].ResType, "DC")==0 && strcmp(atoms_other_chain[k-i].ResType, "DA")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DC")!=0 && strcmp(atoms_other_chain[k-i].ResType, "DT")!=0) ) == FALSE)
			{
				printf("Error\nAtom %s:%s (%u) isn't complement to atom %s:%s (%u). DNA is defect!\n", atoms1_P[i].ResNumber, atoms1_P[i].ResType, i, atoms_other_chain[k-i].ResNumber, atoms_other_chain[k-i].ResType, k-i);
				//printf(max_score, "Error\nDNA is defect! Can not find complement nucleotide.\n");
				exit(1);
			}
		}
	}

	(*compl) = (n_P1+k-1);
	printf("k=%u; compl=%u; n_P1=%u\n", k, (*compl), n_P1 );
}


unsigned int run_find_compl(struct atom *atoms1_P, unsigned int n_P, unsigned int *compl, unsigned int *first_chain_length, int *compl_pair)
{
	unsigned i, k, c, count, k_min, n_P1, n_P2;
	char res1[5], res2[5];

	n_P1 = 0; n_P2 = 0;
	for (i=1; i<=n_P; i++)
	{
		//printf("%c.%s\n", atoms1_P[i].Chain, atoms1_P[i].ResNumber);
		if (atoms1_P[i].Chain == atoms1_P[1].Chain)
		{n_P1 += 1; }
		else {n_P2 += 1; }
	}

	(*first_chain_length) = n_P1;

	sprintf(res1, "%d", compl_pair[1]);
	sprintf(res2, "%d", compl_pair[2]);
	//printf("res1=%s res2=%s\n", res1, res2);

	for (k=2; k<=n_P; k++)
		for (i=1; i<=n_P1; i++)
		{
			if ( (k-i>=1) && (k-i<=n_P2) )
			{
				//printf("k=%u i=%i %s=%s %s=%s\n", k, i, atoms1_P[i].ResNumber, res1, atoms1_P[n_P1+k-i].ResNumber, res2);
				if ( (strcmp(atoms1_P[i].ResNumber, res1)==0) && (strcmp(atoms1_P[n_P1+k-i].ResNumber, res2)==0) )
				{
					//printf("n_P1=%u, k=%u\n", n_P1, k);
					(*compl) = (n_P1+k-1);
					return 0;
				}
			}
		}
	return 1;
}


void BestDiag(double **measures, unsigned int n, unsigned int m, double *S_max,
		unsigned int *i_max, unsigned int *j_max, unsigned int *i_start,
		unsigned int *j_start, unsigned int *i_max_measure, unsigned int *j_max_measure,
		struct atom *atoms1, unsigned int *list_P1, struct atom *atoms2, unsigned int *list_P2,
		unsigned int compl1, unsigned int compl2, unsigned int n_1chain, unsigned int m_1chain){
    //n for n_P1; m for n_P2
    //i_max, j_max for (i,j)-pair with max S-sum
    //i_max_measure, j_max_measure for (i,j)-pair with max measure
    int d; //diagonal number: [-n_P1+1; max(n_P2-n_P1, 0)]
    unsigned int i,j,c, max_i, max_j, max_measure_i, max_measure_j;
    double **S; //S-table
    double measure_max=-1, max_sum=-1;

    //S_max search
    S = (double **)malloc( sizeof(double *)*(n+1) );
    for (i=1; i<=n; i++)  S[i] = (double *)malloc( sizeof(double)*(m+1) );


    //printf("n=%u\n", n);
    for (d=(int)(-n+1); d<=(int)(m_1chain-1); d++){
    	for (i = ((-d+1 > 1) ? -d+1 : 1); i <= (n+d<m_1chain ? n : m_1chain-d); i++){
    	//for i=max(-d+1, 1) to (n_P1 if n_P1+d<m_1chain and m_1chain otherwise)
    		j = d+i;
    		//printf("i=%u j=%u Res1№=%s chain1=%c Res2№=%s chain2=%c\n", i, j, atoms1[list_P1[i]].ResNumber, atoms1[list_P1[i]].Chain, atoms2[list_P2[j]].ResNumber, atoms2[list_P2[j]].Chain);
    		//printf("compl1-i+1=%u compl2-j+1=%u Res1№=%s chain1=%c Res2№=%s chain2=%c\n", compl1-i+1, compl2-j+1, atoms1[list_P1[compl1-i+1]].ResNumber, atoms1[list_P1[compl1-i+1]].Chain, atoms2[list_P2[compl2-j+1]].ResNumber, atoms2[list_P2[compl2-j+1]].Chain);

    		S[i][j] = measures[i][j] + ((i==1 || j==1 || i==n_1chain+1 ) ? 0 : S[i-1][j-1]);
    		//if ( (compl1-i+1)>0 && (compl2-j+1)>0 && (compl1-i+1)<=n_1chain*2 && (compl2-j+1)<=m_1chain*2 )
    		//	{S[i][j] += measures[compl1-i+1][compl2-j+1];}
    		if ( (compl1-i+1)>0 && (compl2-j+1)>0 && (compl1-i+1)<=n && (compl2-j+1)<=m )
    			{S[i][j] += measures[compl1-i+1][compl2-j+1];}
    		if (S[i][j] < 0)  S[i][j] = 0;
    		if (S[i][j] > max_sum){
    			max_sum = S[i][j];
    			max_i=i; max_j=j;
    		}
    		c += 1;
    	}
    }
    (*S_max) = max_sum;
    (*i_max) = max_i;
    (*j_max) = max_j;


    //reverse move
    i=max_i; j=max_j;
    measure_max = 0;
    max_measure_i = i-2;
    max_measure_j = j+2;
    while (i>0 && j>0 && S[i][j]>0){
    	if (measures[i][j] > measure_max){
    		measure_max = measures[i][j];
    		max_measure_i = i;
    		max_measure_j = j;
    	}
    	if ( (compl1-i+1)>0 && (compl2-j+1)>0 && (compl1-i+1)<=n && (compl2-j+1)<=m && measures[compl1-i+1][compl2-j+1] > measure_max ){
    		measure_max = measures[compl1-i+1][compl2-j+1];
    		max_measure_i = compl1-i+1;
    		max_measure_j = compl2-j+1;
    	}
    	i--;
    	j--;
    }

    (*i_max_measure) = max_measure_i;
    (*j_max_measure) = max_measure_j;

    *i_start = i+1;
    *j_start = j+1;

    //Enable in test mode

    //print S-table
    /*for (j=1; j<=m_1chain; j++) printf("%4d", j);
    puts("");
    for (i=n; i>=1; i--){
    	printf("%2d", i);
    	for (j=1; j<=m_1chain; j++){
    		printf("%3.0f ", S[i][j]);
    	}
    	puts("");
    }
    for (j=1; j<=m_1chain; j++) printf("%4d", j);
    puts("");*/
}

/***************************
*		WRITE TO PDB
****************************/

unsigned int createPDB(char *filename, char *title){
  FILE * flow_out;

  flow_out = fopen(filename,"w");
  if (flow_out == NULL) {
    perror(filename);
    exit(1);
  }
  fprintf(flow_out,"HEADER%60s\n",title);
  fprintf(flow_out,"TITLE%15s\n",title);
  fclose(flow_out);

  return 0;
  }

unsigned int writetoPDB(char *filename,
                struct atom *atoms, unsigned int natoms){
  FILE * flow_out;
  unsigned int i;

  flow_out = fopen(filename, "a");
  if (flow_out == NULL) {
    perror(filename);
    exit(1);
  }

  i=1;
  //If first line in output.pdb is xxx
  //printf("ATOM %5s %8.3f%8.3f%8.3f\n", atoms[i].AtomNumber, atoms[i].XCoord, atoms[i].YCoord, atoms[i].ZCoord);

  for (i=1; i<=natoms; i++){
    fprintf(flow_out, "ATOM  %5s  %-3s %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f          %c  \n",
atoms[i].AtomNumber,
atoms[i].AtomName, atoms[i].ResType, atoms[i].ChainRenamed, atoms[i].ResNumber,
atoms[i].XCoord, atoms[i].YCoord, atoms[i].ZCoord, atoms[i].Occup, atoms[i].B,
atoms[i].AtomName[0]);
}
fclose(flow_out);
return 0;
}

unsigned int copyPDB(char *filename, char* src){
  FILE * flow_out;
  flow_out = fopen(filename, "a");
  if (flow_out == NULL) {
    perror(filename);
    exit(1);
  }
  unsigned int m, n, w;
  struct atomname* list;
  struct atom* dna, *prot, *wat;
  char *prot_chains = (char *)malloc( sizeof(char)*(2+1) );
  char *dna_chains = (char *)malloc( sizeof(char)*(2+1) );
  readerPDB(src, &m, &n, &w, 1024, &dna, &dna_chains, &prot, &prot_chains, &wat, &list);
  writetoPDB(filename, dna, m);
  writetoPDB(filename, prot, n);
  fclose(flow_out);
  return 0;
}

unsigned int endPDB(char *filename){
  FILE * flow_out;

  flow_out = fopen(filename,"a");
  if (flow_out == NULL) {
    perror(filename);
    exit(1);
  }

  fprintf(flow_out,"ENDMDL\nEND\n");
  fclose(flow_out);
  return 0;
  }
/***************************
*		end	WRITE TO PDB
****************************/

/*****************************************
*            distances() print mode
*             instead of write mode (check further)
******************************************/
unsigned int distances_print ( double *** dlist, char *atomtype1, char *moltype1,
                         char *atomtype2, char *moltype2,
                         unsigned int dnanum, unsigned int protnum, unsigned int waternum,
                         struct atom *atoms_dna, struct atom *atoms_prot, struct atom *atoms_wat,
                         unsigned int *a1_pros, unsigned int *a2_pros, char *outfile)
{
  FILE * flow_out;
  unsigned int i, j, k, num1, ii;
  struct geomvector cont;

  flow_out = fopen(outfile,"w");
  if (flow_out == NULL) {
    perror(outfile);
    exit(1);
  }

  fprintf(flow_out, "Nuc.\tNo.\tCh.\tAtom\tX\tY\tZ\tDistances to %s atoms of %s molecule\n",atomtype2,(strcmp(moltype2,"p")==0)?"protein":"dna");

  if(strcmp(moltype1,"p") == 0) num1 = protnum;
  else if(strcmp(moltype1,"d") == 0) num1 = dnanum;

  (*a1_pros)=0;

  if( (strcmp(moltype1,"p") == 0) && (strcmp(moltype2,"p") == 0) )
    {
	(*dlist) = (double **) malloc (protnum * sizeof(double));
    for (j = 0; j < protnum; j++)
      (*dlist)[j] = (double *) malloc (protnum * sizeof(double));
  }
  if( (strcmp(moltype1,"p") == 0) && (strcmp(moltype2,"d") == 0) )
    {
	(*dlist) = (double **) malloc (protnum * sizeof(double));
    for (j = 0; j < protnum; j++)
      (*dlist)[j] = (double *) malloc (dnanum * sizeof(double));
  }
  if( (strcmp(moltype1,"d") == 0) && (strcmp(moltype2,"p") == 0) )
    {
	(*dlist) = (double **) malloc (dnanum * sizeof(double));
    for (j = 0; j < dnanum; j++)
      (*dlist)[j] = (double *) malloc (protnum * sizeof(double));
  }
  if( (strcmp(moltype1,"d") == 0) && (strcmp(moltype2,"d") == 0) )
    {
	(*dlist) = (double **) malloc (dnanum * sizeof(double));
    for (j = 0; j < dnanum; j++)
      (*dlist)[j] = (double *) malloc (dnanum * sizeof(double));
  }

  if ((*dlist) == NULL)
  {
    perror("malloc");
    exit(1);
  }

  for (i=1; i<num1; i++)
  {
    if(strcmp(moltype1,"p") == 0)
    {
      if(strcmp(atoms_prot[i].AtomName, atomtype1) == 0)
      {
        fprintf(flow_out,"%s\t%s%c\t%c\t%s\t",
          atoms_prot[i].ResType, atoms_prot[i].ResNumber, atoms_prot[i].Insertion,
          atoms_prot[i].Chain, atoms_prot[i].AtomName);
        fprintf(flow_out, "%f\t", atoms_prot[i].XCoord);
        fprintf(flow_out, "%f\t", atoms_prot[i].YCoord);
        fprintf(flow_out, "%f\t", atoms_prot[i].ZCoord);

        if(strcmp(moltype2,"p") == 0)
        {
          (*a2_pros)=0;
              //printf("%u ",*a1_pros); // Enable in test mode
          for (j=1; j<=protnum; j++)
          {
            if(strcmp(atoms_prot[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_prot[i]), location(atoms_prot[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }
          }
        else if(strcmp(moltype2,"d") == 0)
        {
         (*a2_pros)=0;
          for (j=1; j<=dnanum; j++)
          {
            if(strcmp(atoms_dna[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_prot[i]), location(atoms_dna[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }  //end for
          }  // end if
        fprintf(flow_out, "\n");
	    (*a1_pros)++;
        }  // end if
      }  // end if

    else if(strcmp(moltype1,"d") == 0)
    {
      if(strcmp(atoms_dna[i].AtomName, atomtype1) == 0)
      {
        fprintf(flow_out,"%s\t%s%c\t%c\t%s\t",
          atoms_dna[i].ResType, atoms_dna[i].ResNumber, atoms_dna[i].Insertion,
          atoms_dna[i].Chain, atoms_dna[i].AtomName);
        fprintf(flow_out, "%f\t", atoms_dna[i].XCoord);
        fprintf(flow_out, "%f\t", atoms_dna[i].YCoord);
        fprintf(flow_out, "%f\t", atoms_dna[i].ZCoord);

        if(strcmp(moltype2,"p") == 0)
        {
         (*a2_pros)=0;
          for (j=1; j<=protnum; j++)
          {
            if(strcmp(atoms_prot[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_dna[i]), location(atoms_prot[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }  // end for
          }  //  end if
        else if(strcmp(moltype2,"d") == 0)
        {
          (*a2_pros)=0;
          for (j=1; j<=dnanum; j++)
          {
            if(strcmp(atoms_dna[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_dna[i]), location(atoms_dna[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }  // end for
          }  //end if
        fprintf(flow_out, "\n");
        (*a1_pros)++;
        }  // end if
      }  // end if
    }  // end for
  fclose(flow_out);
  return 0;
  }

/*******************
* distances():
*  compute distances between all atoms of given types (moltype = 'p' or 'd')
*
*******************/


unsigned int distances ( double *** dlist, char *atomtype1, char *moltype1,
                         char *atomtype2, char *moltype2,
                         unsigned int dnanum, unsigned int protnum, unsigned int waternum,
                         struct atom *atoms_dna, struct atom *atoms_prot, struct atom *atoms_wat, //char *outfile,
                         unsigned int *a1_pros, unsigned int *a2_pros)
{
  //FILE * flow_out;
  unsigned int i, j, k, num1, ii;
  struct geomvector cont;

  if(strcmp(moltype1,"p") == 0) num1 = protnum;
  else if(strcmp(moltype1,"d") == 0) num1 = dnanum;

  (*a1_pros)=0;

  if( (strcmp(moltype1,"p") == 0) && (strcmp(moltype2,"p") == 0) )
    {
	(*dlist) = (double **) malloc (((protnum+1)*(protnum+1)+1) * sizeof(double));
    for (j = 0; j <= protnum; j++)
      (*dlist)[j] = (double *) malloc ((protnum+1) * sizeof(double));
  }
  if( (strcmp(moltype1,"p") == 0) && (strcmp(moltype2,"d") == 0) )
    {
	(*dlist) = (double **) malloc (((dnanum+1)*(protnum+1)+1) * sizeof(double));
    for (j = 0; j <= protnum; j++)
      (*dlist)[j] = (double *) malloc ((dnanum+1) * sizeof(double));
  }
  if( (strcmp(moltype1,"d") == 0) && (strcmp(moltype2,"p") == 0) )
    {
	(*dlist) = (double **) malloc (((protnum+1)*(dnanum+1)+1) * sizeof(double));
    for (j = 0; j <= dnanum; j++)
      (*dlist)[j] = (double *) malloc ((protnum+1) * sizeof(double));
  }
  if( (strcmp(moltype1,"d") == 0) && (strcmp(moltype2,"d") == 0) )
    {
	(*dlist) = (double **) malloc (((dnanum+1)*(dnanum+1)+1) * sizeof(double));
    for (j = 0; j <= dnanum; j++)
      (*dlist)[j] = (double *) malloc ((dnanum+1) * sizeof(double));
  }

  if ((*dlist) == NULL)
  {
    perror("malloc");
    exit(1);
  }

  for (i=1; i<=num1; i++)
  {
    if(strcmp(moltype1,"p") == 0)
    {
      if(strcmp(atoms_prot[i].AtomName, atomtype1) == 0)
      {
        /*fprintf(flow_out,"%s\t%s%c\t%c\t%s\t",
          atoms_prot[i].ResType, atoms_prot[i].ResNumber, atoms_prot[i].Insertion,
          atoms_prot[i].Chain, atoms_prot[i].AtomName);
        fprintf(flow_out, "%f\t", atoms_prot[i].XCoord);
        fprintf(flow_out, "%f\t", atoms_prot[i].YCoord);
        fprintf(flow_out, "%f\t", atoms_prot[i].ZCoord);*/

        if(strcmp(moltype2,"p") == 0)
        {
          (*a2_pros)=0;
              //printf("%u ",*a1_pros);
          for (j=1; j<=protnum; j++)
          {
            if(strcmp(atoms_prot[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_prot[i]), location(atoms_prot[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              //fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }
          }
        else if(strcmp(moltype2,"d") == 0)
        {
         (*a2_pros)=0;
          for (j=1; j<=dnanum; j++)
          {
            if(strcmp(atoms_dna[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_prot[i]), location(atoms_dna[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              //fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }  //end for
          }  // end if
        //fprintf(flow_out, "\n");
	    (*a1_pros)++;
        }  // end if
      }  // end if

    else if(strcmp(moltype1,"d") == 0)
    {
      if(strcmp(atoms_dna[i].AtomName, atomtype1) == 0)
      {
        /*fprintf(flow_out,"%s\t%s%c\t%c\t%s\t",
          atoms_dna[i].ResType, atoms_dna[i].ResNumber, atoms_dna[i].Insertion,
          atoms_dna[i].Chain, atoms_dna[i].AtomName);
        fprintf(flow_out, "%f\t", atoms_dna[i].XCoord);
F        fprintf(flow_out, "%f\t", atoms_dna[i].YCoord);
        fprintf(flow_out, "%f\t", atoms_dna[i].ZCoord);*/

        if(strcmp(moltype2,"p") == 0)
        {
         (*a2_pros)=0;
          for (j=1; j<=protnum; j++)
          {
            if(strcmp(atoms_prot[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_dna[i]), location(atoms_prot[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              //fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }  // end for
          }  //  end if
        else if(strcmp(moltype2,"d") == 0)
        {
          (*a2_pros)=0;
          for (j=1; j<=dnanum; j++)
          {
            if(strcmp(atoms_dna[j].AtomName, atomtype2) == 0)
            {
              cont = fromto(location(atoms_dna[i]), location(atoms_dna[j]));
              (*dlist)[(*a1_pros)][(*a2_pros)] = sqrt(scalar(cont,cont));
              //fprintf(flow_out, "%f\t", (*dlist)[(*a1_pros)][(*a2_pros)]);
			  (*a2_pros)++;
              }
            }  // end for
          }  //end if
        //fprintf(flow_out, "\n");
        (*a1_pros)++;
        }  // end if
      }  // end if
    }  // end for
  //fclose(flow_out);
  return 0;
  }

/******************************************
*            readerPDB
*******************************************/
unsigned int readerPDB (char *filename, unsigned int *dnanum, unsigned int *protnum,
              unsigned int *watnum,
              unsigned int maxnumber, struct atom **patoms_dna, char **chains_dna,
              struct atom **patoms_prot, char **chains_prot, struct atom **patoms_wat,
              struct atomname **plist)
{
  FILE *flow_pdb;
  char c[82]; //for string in pdb-file
  char p[10];
  char trash[57];
  char rtype[7];
  char CurrentModel[5];
  char PDBCode[5];
  unsigned int i, j, k;
  int atomflag, changemodel; //atomflag changes if PDB-model changes
  unsigned int maxlist = 128;
  unsigned int maxnumber1 = maxnumber;
  unsigned int maxnumber2 = maxnumber;
  unsigned int maxnumber3 = maxnumber;
  unsigned int recordlength;
  unsigned int dna_chains_num = 0;
  unsigned int dna_chains_max = 2;
  unsigned int prot_chains_num = 0;
  unsigned int prot_chains_max = 2;
  unsigned int result; //for search current atom chain in chain-list
  unsigned int count = 0; //counter

  struct atom currentatom;

  (*dnanum) = 0;
  (*protnum) = 0;
  (*watnum) = 0;

  (*patoms_dna) = (struct atom *) malloc ((maxnumber1+1) * sizeof(struct atom));
  (*patoms_prot) = (struct atom *) malloc ((maxnumber2+1) * sizeof(struct atom));
  (*patoms_wat) = (struct atom *) malloc ((maxnumber3+1) * sizeof(struct atom));

  rtype[6]='\0';
  /* Opening PDB file */
  flow_pdb = fopen(filename,"r");
  if (flow_pdb == NULL)
  {
    perror (filename);
    exit(1);
  }

  currentatom.PDBCode = (char *)malloc( sizeof(char)*6 );
  currentatom.ModelNumber = (char *)malloc( sizeof(char)*6 );
  currentatom.AtomNumber = (char *)malloc( sizeof(char)*7 );
  currentatom.AtomName = (char *)malloc( sizeof(char)*6 );
  currentatom.ResType = (char *)malloc( sizeof(char)*4 );
  currentatom.ResNumber = (char *)malloc( sizeof(char)*7 );

  strcpy(currentatom.PDBCode,"0XXX");
  strcpy(currentatom.ModelNumber,"0");
  strcpy(CurrentModel,"0");


  fgets (c, 80, flow_pdb);
  sscanf(c, "%6c%56c%4c",rtype,trash, PDBCode);
  PDBCode[4] = 0;
  if(strcmp(rtype,"HEADER") == 0) {
    strcpy(currentatom.PDBCode, PDBCode);
  }

  atomflag = TRUE;
  while( atomflag && (strcmp(currentatom.ModelNumber, CurrentModel) == 0 ) ){
    atomflag = readatom(flow_pdb, &currentatom, &changemodel, 1);

    	if( changemodel ) {
      		if (strcmp(currentatom.ModelNumber,"1") == 0) {
        		strcpy(CurrentModel, "1");
      		}
    	}


	    if( strcmp(currentatom.ModelNumber, CurrentModel) == 0 ) {
	    //printf("%s, %s\n", currentatom.ResType, currentatom.ResNumber);
		  if ( ifprot(currentatom.ResType) == TRUE ) {
		     (*protnum)++;
		     if( (*protnum) >= maxnumber2 ){
		        maxnumber2 += maxnumber;
		        (*patoms_prot) = (struct atom *) realloc ((*patoms_prot), maxnumber2 * sizeof(struct atom));
		     }
		     count = 0;
		     result = 0;
		     while (++count != prot_chains_num+1){
		     	if ((*chains_prot)[count] == currentatom.Chain)
		        	result = 1;
		     }
		     if (result == 0){
		     	(*chains_prot)[count] = currentatom.Chain;
		     	prot_chains_num++;
		     	if (prot_chains_num >= prot_chains_max){
		     		prot_chains_max++;
		     		(*chains_prot) = (char *)realloc( (*chains_prot), sizeof(char)*(prot_chains_max+1) );
		     	}
		     	(*chains_prot)[prot_chains_num+1] = '\0';
		     }
		     atomcpy ( &((*patoms_prot)[(*protnum)]), currentatom);
		  }
		  else if ( ifdna(currentatom.ResType) == TRUE ) {

		     (*dnanum)++;
		     if( (*dnanum) >= maxnumber1 ){
			maxnumber1 += maxnumber;
			(*patoms_dna) = (struct atom *) realloc ((*patoms_dna), maxnumber1 * sizeof(struct atom));
		     }
		     count = 0;
		     result = 0;
		     while (++count != dna_chains_num+1){
		     	if ((*chains_dna)[count] == currentatom.Chain)
				result = 1;
		     }
		     if (result == 0){
			     printf("ADD DNA CHAIN %c for %s\n", currentatom.Chain, currentatom.PDBCode);
		     	(*chains_dna)[count] = currentatom.Chain;
		     	dna_chains_num++;
		     	if (dna_chains_num >= dna_chains_max){
		     		dna_chains_max += 2;
		     		(*chains_dna) = (char *)realloc( (*chains_dna), sizeof(char)*(dna_chains_max+1) );
		     	}
		     	(*chains_dna)[dna_chains_num+1] = '\0';
		     }
		     atomcpy ( &((*patoms_dna)[(*dnanum)]), currentatom);
		     if ( (*dnanum) > 1 )
		     if ( strcmp(currentatom.AtomName, (*patoms_dna)[(*dnanum)-1].AtomName ) == 0 && strcmp(currentatom.ResNumber, (*patoms_dna)[(*dnanum)-1].ResNumber) == 0 )
    		     	{

    		     		(*dnanum)--;
		        } // if AtomName and ResNumber differs from previous DNA atom to avoid atoms with alternative position
		  }
		  else if ( strcmp(currentatom.ResType,"WAT") == 0 ||
		            strcmp(currentatom.ResType,"HOH") == 0 ||
		            strcmp(currentatom.ResType,"H2O") == 0 ) {
		     (*watnum)++;
		     if( (*watnum) >= maxnumber3 ){
		        maxnumber3 += maxnumber;
		        (*patoms_wat) = (struct atom *) realloc ((*patoms_wat), maxnumber3 * sizeof(struct atom));
		     }
		     atomcpy ( &((*patoms_wat)[(*watnum)]), currentatom);
		  } /* if the atom is from protein / else if it is from dna / else it is water */
	   } /* if first model or no models */

  } /* while */

  fclose (flow_pdb);    /*close infile*/

  return (*dnanum + *protnum);
} /* readerPDB */


/***************************************************************
*          ifprot
* Returns TRUE if the string is a name of an aminoacid residue;
* returns FALSE otherwise.
****************************************************************/

unsigned int ifprot(char *res)
{
  unsigned int result = FALSE;
  if ( strcmp(res, "ALA") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "ARG") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "ASN") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "ASP") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "CYS") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "GLN") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "GLU") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "GLY") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "HIS") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "ILE") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "LEU") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "LYS") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "MET") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "PHE") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "PRO") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "SER") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "THR") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "TRP") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "TYR") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "VAL") == 0 ) {
    result = TRUE;
  }

  return result;
} /* ifprot */

/******************************************************************
*          ifdna
* Returns TRUE if the string is a name of a nucleotide residue;
* returns FALSE otherwise.
*******************************************************************/
unsigned int ifdna(char *res)
{
  unsigned int result = FALSE;
  char nt;
  if ( strcmp(res, "A") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "T") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "G") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "C") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "U") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "DA") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "DT") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "DG") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "DC") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "5IT") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "+A") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "+C") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "+G") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "+U") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "PST") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "3DR") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "ASU") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "ORP") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "IMP") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "NRI") == 0 ) {
    result = TRUE;
  }
  if ( strcmp(res, "YRR") == 0 ) {
    result = TRUE;
  }

  if (result == FALSE) {
    result = 1 - convertNT(res, &nt);
  }

  return result;
} /* ifdna */

/*****************************************************
* Readatom: reads the next atom from the file;
* if key = 0 - any atom, if key = 1 - any heavy atom,
* key = 2 - protein heavy atom, key = 3 - DNA/RNA heavy atom
*******************************************************/
unsigned int readatom (FILE *flow_pdb, struct atom *current,
                       int *changemodel, int key)
{
  char c[80];
  char rtype[7], cnatom[6], blank[7],
       atnm[5], alter[2], resnm[4],
       chnid[2], cnres[5], insert[2],
       cx[9], cy[9], cz[9], coccup[7], ctf[7];
  int atomflag;
  unsigned int i,j;

  rtype[6]='\0';
  cnatom[5]='\0';
  atnm[4]='\0';
  alter[1]='\0';
  resnm[3]='\0';
  chnid[1]='\0';
  cnres[4]='\0';
  insert[1]='\0';
  cx[8]='\0';
  cy[8]='\0';
  cz[8]='\0';
  coccup[6]='\0';
  ctf[6]='\0';

  *changemodel = FALSE;
  atomflag = FALSE;

  while ( (atomflag == FALSE) && !feof(flow_pdb) ){

    fgets (c, 80, flow_pdb);
    sscanf(c, "%6c", rtype);

    if(strcmp(rtype,"MODEL ") == 0) {
       *changemodel = TRUE;
       j = 0;
       for (i = 11;i<=14;i++) {
         if ( (int)c[i] > 47 ) {
           current -> ModelNumber[j] = c[i];
           j++;
         }
       }
       current -> ModelNumber[j] = '\0';
    }

    if(strcmp(rtype,"ATOM  ") == 0  || strcmp(rtype,"HETATM") == 0)
    {
      sscanf(c, "%6c%5c%c%4c%c%3c%c%c%4c%c%3c%8c%8c%8c%6c%6c",
      		rtype, cnatom, blank, atnm, alter, resnm, blank, chnid,
      		cnres, insert, blank, cx, cy, cz, coccup, ctf);

      sscanf(atnm, "%s", atnm);
      if (key == 0 ||
          (atnm[0] != 'H' && atnm[0] != '1' &&
           atnm[0] != '2' && atnm[0] != '3') )
      { /* not hydrogen or we accept hydrogens */
        sscanf(resnm, "%s", resnm);
        if( ( (ifprot(resnm) == TRUE) && (key == 2) ) ||
            ( (ifdna(resnm) == TRUE) && (key == 3) ) ||
            ( key == 1 ) || (key == 0) ){
          sscanf(cnatom, "%s", cnatom);
          sscanf(cnres, "%s", cnres);

/*        if (alter[0]==' ') alter[0]='~';
          else printf ("\nAlter code %c in atom %s (%s%s.%s)!",
                        alter[0],cnatom,resnm,cnres,atnm); */


          atomflag = TRUE;
          strcpy(current -> AtomNumber, cnatom);
          strcpy(current -> AtomName, atnm);
          current -> Alter = alter[0];
          strcpy(current -> ResType, resnm);
          if (chnid[0]==' ') chnid[0]='*';
          current -> Chain = chnid[0];
          strcpy(current -> ResNumber, cnres);
          current -> Insertion = insert[0];
          current -> XCoord = atof(cx);
          current -> YCoord = atof(cy);
          current -> ZCoord = atof(cz);
          current -> Occup = atof(coccup);
          current -> B = atof(ctf);
        } /* if the residue name is in accordance with the key */
      } /* if not hydrogen or key is 0 */
    } /* if ATOM || HETATM */
  } /* while */

  return atomflag;
} /* readatom */


/**********************
* atomcpy():
* Copy all properties of given atom to the target atom
**********************/

void atomcpy(struct atom *target, struct atom source)
{
  target -> PDBCode =
    (char *)malloc( sizeof(char)*(strlen(source.PDBCode)+1) );
  target -> ModelNumber =
    (char *)malloc( sizeof(char)*(strlen(source.ModelNumber)+1) );
  target -> AtomNumber =
    (char *)malloc( sizeof(char)*(strlen(source.AtomNumber)+1) );
  target -> AtomName =
    (char *)malloc( sizeof(char)*(strlen(source.AtomName)+1) );
  target -> ResType =
    (char *)malloc( sizeof(char)*(strlen(source.ResType)+1) );
  target -> ResNumber =
    (char *)malloc( sizeof(char)*(strlen(source.ResNumber)+1) );

  strcpy(target -> PDBCode, source.PDBCode);
  strcpy(target -> ModelNumber, source.ModelNumber);
  strcpy(target -> AtomNumber, source.AtomNumber);
  strcpy(target -> AtomName, source.AtomName);
  strcpy(target -> ResType, source.ResType);
  strcpy(target -> ResNumber, source.ResNumber);
  target -> Alter = source.Alter;
  target -> Chain = source.Chain;
  target -> Insertion = source.Insertion;
  target -> XCoord = source.XCoord;
  target -> YCoord = source.YCoord;
  target -> ZCoord = source.ZCoord;
  target -> Occup = source.Occup;
  target -> B = source.B;
  return;
}

void deleteAtom(struct atom* target) {
	free(target->PDBCode);
	free(target->ModelNumber);
	free(target->AtomNumber);
	free(target->AtomName);
	free(target->ResType);
	free(target->ResNumber);
}

/**********************
* atomlistcpy():
* Copy all atoms of given list of atoms to the target list
**********************/

void atomlistcpy(struct atom **target, struct atom *source, unsigned int leng)
{
	unsigned int i;
	(*target) = (struct atom *)malloc( sizeof(struct atom)*(leng+1) );
	for (i=1; i<=leng; i++)
		atomcpy(&(*target)[i], source[i]);
}


void atomlistmerge(struct atom **target, unsigned int *len, struct atom *source1, unsigned int len1, struct atom *source2, unsigned int len2)
{
	unsigned int i, c;
	(*target) = (struct atom *)malloc( sizeof(struct atom)*(len1+len2+1) );
	(*len) = len1 + len2;
	for (i=1; i<=len1; i++)
		atomcpy(&(*target)[i], source1[i]);

	c = len1+1;
	for (i=1; i<=len2; i++, c++)
		atomcpy(&(*target)[c], source2[i]);
}

void deleteAtomlist(struct atom* target, unsigned int len) {
	for (int i=1;i<=len;i++) {
		deleteAtom(&target[i]);
	}
	free(target);
}
/************************
* location
************************/
struct geompoint location(struct atom a)
{
  struct geompoint result;
  result.X = a.XCoord;
  result.Y = a.YCoord;
  result.Z = a.ZCoord;
  return result;
} /* location */

#define EPSILON 0.000001

/**********************************
*  Scalar product of two vectors
**********************************/
double scalar (struct geomvector v, struct geomvector u)
{
  return v.X*u.X + v.Y*u.Y + v.Z*u.Z;
} /* scalar */

/***********************************************
* "fromto" - the vector from point A to point B
************************************************/
struct geomvector fromto(struct geompoint A, struct geompoint B)
{
  struct geomvector result;
  result.X = B.X - A.X;
  result.Y = B.Y - A.Y;
  result.Z = B.Z - A.Z;
  return result;
} /* geomvector */

/**********************************************************
* Square distance between point A and straight line (B,V)
***********************************************************/
double sqdistance(struct geompoint A, struct geompoint B,
				  struct geomvector V)
{
  struct geomvector ab;
  double sqab, vab, sqv, result;

  ab = fromto(A,B);

  vab = scalar(V,ab);
  sqab = scalar(ab,ab);
  sqv = scalar(V,V);

  if (sqv > 0)
  {
    result = sqab - vab*vab/sqv;
  }
  else { /* V has zero length */
    perror("Straight line with zero direction!!!");
    exit(1);
  }
  return result;
} /* sqdistance */

/******************************************
* "num_x_vec" - a vector V multiplied by
*  a real number "coef"
******************************************/
struct geomvector num_x_vec(double coef, struct geomvector V)
{
  struct geomvector result;
  result.X = coef*V.X;
  result.Y = coef*V.Y;
  result.Z = coef*V.Z;
  return result;
} /* num_x_vec */

/******************************************
*  "pointplusvec" - the shift of a point A
*  with a vector V
*******************************************/
struct geompoint pointplusvec(struct geompoint A, struct geomvector V)
{
  struct geompoint result;
  result.X = A.X + V.X;
  result.Y = A.Y + V.Y;
  result.Z = A.Z + V.Z;
  return result;
} /* pointplusvec */

/************************
* the sum of two vectors
*************************/
struct geomvector sumvec(struct geomvector V, struct geomvector W)
{
  struct geomvector result;
  result.X = V.X + W.X;
  result.Y = V.Y + W.Y;
  result.Z = V.Z + W.Z;
  return result;
} /* sumvec */

/*****************************************
* "nearest" - the point on the line (A,V)
* that is nearest to the line (B,W)
******************************************/
struct geompoint nearest(struct geompoint A, struct geomvector V,
                         struct geompoint B, struct geomvector W)
{
  struct geompoint result;
  struct geomvector ortV, ortW;
  double x,s,t;

  if( (scalar(V,V) > 0) && (scalar(W,W) > 0)){
    ortV = ort(V);
    ortW = ort(W);
    t = 1 - scalar(ortV,ortW)*scalar(ortV,ortW);
    if ( t*t < EPSILON ){
      perror("Nearest: parallel lines!!!");
      exit(1);
    }
    s = scalar(fromto(A,B),sumvec(ortV,num_x_vec((-1)*scalar(ortV,ortW),ortW)));
    x = s/t;
    result = pointplusvec(A,num_x_vec(x,ortV));
  }
  else { /* V or W has zero length */
    perror("Straight line with zero direction!!!");
    exit(1);
  }
  return result;
} /* nearest */

/*******************
* vector product
*******************/
struct geomvector vecproduct(struct geomvector V, struct geomvector W)
{
  struct geomvector result;
  result.X = V.Y*W.Z - W.Y*V.Z;
  result.Y = V.Z*W.X - W.Z*V.X;
  result.Z = V.X*W.Y - W.X*V.Y;
  return result;
} /* vecproduct */

/***************************************
* Projection of point B onto line (A,V)
****************************************/
struct geompoint proj(struct geompoint B,
                      struct geompoint A, struct geomvector V)
{
  struct geomvector ort_V;
  struct geompoint result;

  ort_V = ort(V);
  result = pointplusvec(A,num_x_vec(scalar(fromto(A,B),ort_V), ort_V));

  return result;
} /* proj */

/************************************
* Derived point = t*A + (1-t)*B
*************************************/
struct geompoint derpoint(double t, struct geompoint A, struct geompoint B)
{
  struct geompoint result;
  struct geomvector connect;

  connect = fromto(A,B);
  result = pointplusvec(A, num_x_vec(t, connect));

  return result;
} /* derpoint */

/**********************************************
* Vector of unit length in the given direction
***********************************************/
struct geomvector ort(struct geomvector V)
{
  struct geomvector result;
  double length;

  length = sqrt(scalar(V,V));
  if(length > 0) result = num_x_vec(1/length, V);
  else {
   perror("Vector of zero length!!!");
   exit(1);
  }
  return result;
} /* ort */

/**********************************************************
* Turn of the vector V at the angle "turn" around the axis
***********************************************************/
struct geomvector turnvec(double turn, struct geomvector axis, struct geomvector V)
{
  struct geomvector VX, W, Wturn, result;
  struct coordsystem new;

  VX = vecproduct(axis,V);
  if( scalar(VX,VX) == 0 ) result = V; /* Turn around a parrallel axis is identity */
  else {
    new.ort_Z = ort(axis);
    new.ort_X = ort(VX);
    new.ort_Y = vecproduct(new.ort_Z, new.ort_X);

    W = changesystem(V, new, 1);  /* W.X is 0 by the construction */
    Wturn.X = - sin(turn)*W.Y;
    Wturn.Y = cos(turn)*W.Y;
    Wturn.Z = W.Z;
    result = changesystem(Wturn, new, -1);
  }
  return result;
} /* turnvec */

/*********************************************
* Changing coordinates systems
* into the new coordinate system (flag = 1)
* or from it (flag = -1)
**********************************************/
struct geomvector changesystem(struct geomvector V, struct coordsystem newcs, int flag)
{
  struct geomvector result;

  if(flag == 0){
    perror("Changesystem: flag must be non-zero!!!");
    result = V;
  }

  if(flag > 0)
  {
    result.X= scalar(V, newcs.ort_X);
    result.Y= scalar(V, newcs.ort_Y);
    result.Z= scalar(V, newcs.ort_Z);
  }
  if(flag < 0)
  {
    result = sumvec(num_x_vec(V.X,newcs.ort_X),
                    sumvec(num_x_vec(V.Y,newcs.ort_Y),num_x_vec(V.Z,newcs.ort_Z)));
  }
  return result;
} /* changesystem */

/*
 * GK file generation
 */
static unsigned int nearestAtom(struct atom from, struct atom* chain, const unsigned int* filter, unsigned int size) {
	double distance = 1.0/0.0;
	unsigned int result = 0;
	for (unsigned int i=0; i< size;i++){
		struct geomvector diff = fromto(location(from), location(chain[filter[i+1]]));
		double cur = scalar(diff, diff);
		if (cur < distance) {
			distance = cur;
			result = i;
		}
	}
	return result;
}

static int imax(int a, int b) {
	return a<b?b:a;
}

static void reverse(char* out) {
	int n = strlen(out);
	for (int i=0, j=n-1; i<j; i++, j--) {
		char tmp = out[i];
		out[i] = out[j];
		out[j] = tmp;
	}
}

static void sputChain(char* out, struct atom* chain, unsigned int* P, unsigned int n_P, char (*coder)(const struct atom*)) {
	for (int i=1;i<=n_P;i++) {
		*(out++) = coder(&chain[P[i]]);
	}
	*out = 0;
}

static void sputChainWithGaps(char* out1, char* out2, struct atom* chain1, unsigned int* P1, unsigned int n_P1, struct atom* chain2, unsigned int* P2, unsigned int n_P2, unsigned int* nearestToCentroid, unsigned int* nearestOnCentroid, char (*coder)(const struct atom*), int useMiddleGaps) {
	int m = n_P1+1;
	int n = n_P2+1;
	double* scoreMatrix = (double*)malloc(sizeof(double) * n * m);
	memset(scoreMatrix, 0, sizeof(double)*n*m);
	char* fromMatrix = (char*)malloc(sizeof(char) * n * m);
	memset(fromMatrix, 0, n*m);
	if (useMiddleGaps) {
		for (int i=1; i<n_P1; i++) {
			for (int j=1; j<n_P2; j++) {
				scoreMatrix[i*n+j] = scoreMatrix[i*n+j-1];
				fromMatrix[i*n+j] = 'L';
				double possibleValue = scoreMatrix[i*n+j-n];
				if (scoreMatrix[i*n+j] <= possibleValue) {
					scoreMatrix[i*n+j] = possibleValue;
					fromMatrix[i*n+j] = 'U';
				}
				struct geomvector diff = fromto(location(chain1[P1[i]]), location(chain2[P2[j]]));
				possibleValue = scoreMatrix[i*n+j-n-1] + /*1.0/ (sqrt(scalar(diff, diff))+0.00001) */(nearestToCentroid[i-1] == j-1 && nearestOnCentroid[j-1] == i-1)+0.0001;
				if (scoreMatrix[i*n+j] <= possibleValue) {
					scoreMatrix[i*n+j] = possibleValue;
					fromMatrix[i*n+j] = '+';
				}
			}
		}
	} else {
		for (int i=1; i<n_P1; i++) {
			for (int j=1; j<n_P2; j++) {
				struct geomvector diff = fromto(location(chain1[P1[i]]), location(chain2[P2[j]]));
				double possibleValue = scoreMatrix[i*n+j-n-1] + /*1.0/ (sqrt(scalar(diff, diff))+0.00001) */(nearestToCentroid[i-1] == j-1 && nearestOnCentroid[j-1] == i-1)+0.0001;
				scoreMatrix[i*n+j] = possibleValue;
				fromMatrix[i*n+j] = '+';
			}
		}
	}
	for (int i=1; i<=n_P1; i++) {
		for (int j=1; j<=n_P2; j++) {
			if (i!=n_P1 && j!=n_P2) continue;
			scoreMatrix[i*n+j] = scoreMatrix[i*n+j-1];
			fromMatrix[i*n+j] = 'L';
			double possibleValue = scoreMatrix[i*n+j-n];
			if (scoreMatrix[i*n+j] <= possibleValue) {
				scoreMatrix[i*n+j] = possibleValue;
				fromMatrix[i*n+j] = 'U';
			}
			struct geomvector diff = fromto(location(chain1[P1[i]]), location(chain2[P2[j]]));
			possibleValue = scoreMatrix[i*n+j-n-1] + /*1.0/ (sqrt(scalar(diff, diff))+0.00001) */(nearestToCentroid[i-1] == j-1 && nearestOnCentroid[j-1] == i-1)+0.0001;
			if (scoreMatrix[i*n+j] <= possibleValue) {
				scoreMatrix[i*n+j] = possibleValue;
				fromMatrix[i*n+j] = '+';
			}
		}
	}/*
	for (int i=0;i<m;i++) {
		for (int j=0;j<n;j++) {
			printf("%c|", fromMatrix[i*n+j]);
		}
		printf("\n");
	}
	for (int i=0;i<m;i++) {
		for (int j=0;j<n;j++) {
			printf("%f|", scoreMatrix[i*n+j]);
		}
		printf("\n");
	}*/
	int i = n_P1, j = n_P2;
	char* out1Ptr = out1, *out2Ptr = out2;
	while (i!=0 || j!=0) {
	//	printf("%c", fromMatrix[i*n+j]);
		if (fromMatrix[i*n+j] == '+') {
			*(out1Ptr++) = coder(&chain1[P1[i]]);
			*(out2Ptr++) = coder(&chain2[P2[j]]);
			i--;
			j--;
		} else if(fromMatrix[i*n+j] == 'U' || j==0) {
			*(out1Ptr++) = coder(&chain1[P1[i]]);
			*(out2Ptr++) = '-';
			i--;
		} else {
			*(out1Ptr++) = '-';
			*(out2Ptr++) = coder(&chain2[P2[j]]);
			j--;
		}
	}
	//printf("\n");
	*out1Ptr = 0;
	*out2Ptr = 0;
	reverse(out1);
	reverse(out2);
	free(fromMatrix);
	free(scoreMatrix);
	//puts(out1);
	//puts("\n");
	//puts(out2);
	//puts("\n");
}

static char encodeDNA(const struct atom* v) {
	return v->ResType[1];
}

static char encodeP(const struct atom* v) {
	if (strcmp(v->ResType, "ALA") == 0) {
		return 'A';
	}
	if (strcmp(v->ResType, "ARG") == 0) {
		return 'R';
	}
	if (strcmp(v->ResType, "ASN") == 0) {
		return 'N';
	}
	if (strcmp(v->ResType, "ASP") == 0) {
		return 'D';
	}
	if (strcmp(v->ResType, "CYS") == 0) {
		return 'C';
	}
	if (strcmp(v->ResType, "GLN") == 0) {
		return 'Q';
	}
	if (strcmp(v->ResType, "GLU") == 0) {
		return 'E';
	}
	if (strcmp(v->ResType, "GLY") == 0) {
		return 'G';
	}
	if (strcmp(v->ResType, "HIS") == 0) {
		return 'H';
	}
	if (strcmp(v->ResType, "ILE") == 0) {
		return 'I';
	}
	if (strcmp(v->ResType, "LEU") == 0) {
		return 'L';
	}
	if (strcmp(v->ResType, "LYS") == 0) {
		return 'K';
	}
	if (strcmp(v->ResType, "MET") == 0) {
		return 'M';
	}
	if (strcmp(v->ResType, "PHE") == 0) {
		return 'F';
	}
	if (strcmp(v->ResType, "PRO") == 0) {
		return 'P';
	}
	if (strcmp(v->ResType, "SER") == 0) {
		return 'S';
	}
	if (strcmp(v->ResType, "THR") == 0) {
		return 'T';
	}
	if (strcmp(v->ResType, "TRP") == 0) {
		return 'W';
	}
	if (strcmp(v->ResType, "TYR") == 0) {
		return 'Y';
	}
	if (strcmp(v->ResType, "VAL") == 0) {
		return 'V';
	}
	return '?';
}

void writegk(const char* path, const char* pdbname, char** infiles, size_t num_pdbs, struct atom** chains, unsigned int* n_chains, unsigned int n, unsigned int centroid_id) {
	FILE* out = fopen(path, "w");
	unsigned int** atoms = malloc(sizeof(unsigned int*) * n);
	unsigned int* n_atoms = malloc(sizeof(unsigned int) * n);
	for (int i=0;i<n;i++) {
		getAtomsNumbers(chains[i], n_chains[i], &atoms[i], &n_atoms[i], i%3==0?"CA":"P");
	}
	unsigned int** nearestToCentroid = malloc(sizeof(unsigned int*) * n);
	unsigned int** nearestOnCentroid = malloc(sizeof(unsigned int*) * n);
	char** formatted = malloc(sizeof(char*) * n);
	char** base = malloc(sizeof(char*) * n);
	char** reference = malloc(sizeof(char*) * 3);
	for (int chain = 0; chain<3; chain++) {
		int compare_to = centroid_id*3+chain;
		unsigned int max_len = 0;
		for (unsigned int i=chain; i<n; i+=3)
			max_len += n_atoms[i];
		reference[chain] = malloc(max_len+1);
		sputChain(reference[chain], chains[compare_to], atoms[compare_to], n_atoms[compare_to], chain==0?encodeP:encodeDNA);
		unsigned int* gaps = malloc(sizeof(unsigned int) * (max_len + 1));
		unsigned int* gapsNeeded = malloc(sizeof(unsigned int) * (max_len + 1)*n/3);
		memset(gaps, 0, sizeof(unsigned int) * (max_len + 1));
		for (unsigned int i=chain; i<n; i+=3) {
			nearestToCentroid[i] = malloc(sizeof(unsigned int) * n_atoms[i]);
			nearestOnCentroid[i] = malloc(sizeof(unsigned int) * n_atoms[compare_to]);
			for (int j=0;j<n_atoms[i];j++) {
				nearestToCentroid[i][j] = nearestAtom(chains[i][atoms[i][j+1]], chains[compare_to], atoms[compare_to], n_atoms[compare_to]);
			}
			for (int j=0;j<n_atoms[compare_to];j++) {
				nearestOnCentroid[i][j] = nearestAtom(chains[compare_to][atoms[compare_to][j+1]], chains[i], atoms[i], n_atoms[i]);
			}
			if (i == compare_to) continue;
			char* buf1 = malloc(strlen(reference[chain]) + n_chains[i] + 1);
			char* buf2 = malloc(strlen(reference[chain]) + n_chains[i] + 1);
			sputChainWithGaps(buf1, buf2, chains[i], atoms[i], n_atoms[i], chains[compare_to], atoms[compare_to], n_atoms[compare_to], nearestToCentroid[i], nearestOnCentroid[i], chain==0?encodeP:encodeDNA, chain==0);
			int gap_id = 0;
			const char* c;
			for (c = buf2; *c != 0; c++) {
				unsigned int gap_len = 0;
				while (*(c++) == '-') {
					gap_len++;
				}
				c--;
				gaps[gap_id] = imax(gap_len, gaps[gap_id]);
				gapsNeeded[gap_id + i/3*(max_len+1)] = gap_len;
				gap_id++;
				if (*c == 0) break;
			}
			unsigned int gap_len = 0;
			while (*(c++) == '-') {
				gap_len++;
			}
			gaps[gap_id] = imax(gap_len, gaps[gap_id]);
			gapsNeeded[gap_id + i/3*(max_len+1)] = gap_len;
			formatted[i] = buf1;
			base[i] = buf2;
		}
		char* tmp = malloc(max_len+1);
		char* itmp=tmp;
		int i;
		for (i=0; reference[chain][i] != 0; i++) {
			for (int j=0;j<gaps[i]; j++) *(itmp++) = '-';
			*(itmp++) = reference[chain][i];
		}
		for (int j=0;j<gaps[i]; j++) *(itmp++) = '-';
		*itmp = 0;
		reference[chain] = tmp;
		formatted[compare_to] = tmp;
		size_t seq_len = strlen(tmp);

		for (int i = chain, j=0; i < n; i+=3, j++) {
			if (i == compare_to) continue;
			tmp = malloc(seq_len+1);
			itmp = tmp;
			int gap_id = 0;
			const char* cb, *cf;
			for (cb = base[i], cf = formatted[i]; *cb != 0; cb++, cf++) {
				size_t have = gapsNeeded[gap_id + i/3*(max_len+1)];
				size_t need = gaps[gap_id];
				if (*cb == '-') {
					if (gap_id == 0 && have < need) {
						for (int k=have; k<need; k++) {
							*(itmp++) = '-';
						}
						gapsNeeded[gap_id + i/3*(max_len+1)] = need;
					}
					*(itmp++) = *cf;
				} else {
					for (int k=have; k<need; k++) {
						*(itmp++) = '-';
					}
					*(itmp++) = *cf;
					gap_id++;
				}
			}
			size_t have = gapsNeeded[gap_id + i/3*(max_len+1)];
			size_t need = gaps[gap_id];
			for (int k=have; k<need; k++) {
				*(itmp++) = '-';
			}
			*(itmp++) = 0;
			formatted[i] = tmp;
		}


		free(gaps);
		free(gapsNeeded);
	}
  printf("In order to avoid confusions, chains in output files were renamed:\n");
  for (int i = 0; i<n;i++) {
    if (i%3==0) {
      fprintf(out, "&%.*s\n", (int) strlen(infiles[i/3]) - 4, infiles[i/3]);
    }
    unsigned int compare_to = i % 3 + centroid_id*3;
    if (i%3==0) {
      printf("Chain %c from %s structure was renamed %c \n", (chains[i])[1].Chain, (chains[i])[1].PDBCode, (chains[i])[1].ChainRenamed);
      printf("Chain %c from %s structure was renamed %c \n", (chains[i+1])[1].Chain, (chains[i+1])[1].PDBCode, (chains[i+1])[1].ChainRenamed);
      printf("Chain %c from %s structure was renamed %c \n", (chains[i+2])[1].Chain, (chains[i+2])[1].PDBCode, (chains[i+2])[1].ChainRenamed);
      fprintf(out, "(%s):%c.%s/1\n", formatted[i+1], (chains[i+1])[1].ChainRenamed, i+1 % 3 == 0?"CA":"P");
      fprintf(out, "(%s):%c.%s/1\n", formatted[i+2], (chains[i+2])[1].ChainRenamed, i+2 % 3 == 0?"CA":"P");
      fprintf(out, "(%s):%c.%s/1\n", formatted[i], (chains[i])[1].ChainRenamed, i % 3 == 0?"CA":"P");
    }
  }
  fclose(out);

}
