/*------------------------------------------------------*/
//BIBLIOTECAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iohb.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

#define IMAX 1000
#define ERRO 0.0001

/*-----------------------------------------------------------------*/
void  le_matriz( char *arquivo, int *M, int *N, int *naozeros, int **colptr, int **linhas, double **valores){
	  
	int retorno, nrhs;
	char *tipo; 
	retorno = readHB_info(arquivo, M, N, naozeros, &tipo, &nrhs);
    
	if (retorno == 0){
	        printf("Erro ao ler as informaçõess da matriz!\n");
		exit(-1);
	}

	printf("Linhas: %d \t Colunas: %d \t Não Zeros: %d \n\n", *M, *N, *naozeros);
	
	*valores = (double *) malloc (*naozeros * sizeof(double));
	*linhas  = (int *) malloc (*naozeros * sizeof(int));
	*colptr  = (int *) malloc ((*N+1) * sizeof(int));

	retorno = readHB_mat_double(arquivo, *colptr, *linhas, *valores);

	     
	if (retorno == 0){
        	printf("Erro ao ler os valores da matriz!\n");
        	exit(-1);
    	}
}
/*-----------------------------------------------------------------*/
void escreve_matriz(int id, int M, int N, int naozeros,int nproc, int *colptr, int *linhas, double *valores){
	int i;
	
	// if(id != 4) return;

	// printf("PTR:\n");
	// for(i=0; i<N+1; i++){
	// 	printf("%d ", colptr[i]);
	// }
	// printf("\n");
	
	// printf("LINHAS id:%d:\n",id);
	// for(i=0; i<naozeros; i++){
	// 	printf("%d ", linhas[i]);
	// }
	// printf("\n\n");
	
	// for(i=0; i<naozeros; i++){
	// 	printf("ID:%d i:%f\n", id,valores[i]);
	// }
	// printf("\n\n");
	/**/
}



void escreve_vetor(int id, double *vet, int n){
	int i;
	printf("ID=%d\n",id);
	for(i=0; i<n; i++){
		printf("% 5.2lf\n", vet[i]);
	}
}

void escreve_vetorInt(int id, int *vet, int n){
	int i;
	printf("ID=%d\n",id);
	for(i=0; i<n; i++){
		printf("%d\n", vet[i]);
	}
}

void inicializa_vetor(int n, double *vet, int valor){
	for(int i=0; i<n; i++){
		vet[i] = valor;
	}
}

void inicializa_vetorInt(int n, int *vet, int valor){
	for(int i=0; i<n; i++){
		vet[i] = valor;
	}
}

void calcula_transposta(int id,double *mat,int *linhas, int *colptr, double *mat_t, int *colunas, int *linptr, int N, int naozeros){
	int i,j,temp,lin,dest,last,cumsum=0;
	inicializa_vetorInt(N+1, linptr,0);
	for(i=0;i<naozeros;i++){
		linptr[linhas[i]-1]++;  
		// linhas[i]-1
	}

	for(i = 0;i<N;i++){
		temp = linptr[i];
		linptr[i] = cumsum;
		cumsum+=temp;
	}
	linptr[N] = naozeros;
	// talvez precise aumentar um número nos valores


	for(i=0;i<N;i++){
		for(j=colptr[i]-1;j<colptr[i+1]-1;j++){
			lin = linhas[j]-1;
			dest = linptr[lin];

			colunas[dest] = i+1;
			mat_t[dest] = mat[j];

			linptr[lin]++;
		}
	}

	last = 0;
	for(i =0;i<=N;i++){
		temp = linptr[i];
		linptr[i] = last + 1;
		last = temp;
	}
}


void produtoMatrizVetor(int n, double *mat, double *vet, double *res){
	int i,j;
	for(i=0;i<n;i++){
		res[i] = 0;
		for(j=0;j<n;j++)
			res[i] += mat[i*n+j] * vet[j];
		
	}	
}

void produto_matriz_vetor (int id, int nproc, double *val, int *lin, int *ptr, double *vet, double *res, int N){
	int i,j,k = 0;
	printf("---produto_matriz_vetor---\n");
	
	printf("ID:%d\n",id);

	if(id != nproc-1){  //id=0;N=5;nproc=5
		for(i=(id*(N/nproc));i<(id+1)*(N/nproc);i++){
			// inic=0;fim=2 
			res[k] = 0;
			for(j=ptr[i]-1;j<ptr[i+1]-1;j++){
				res[k] += val[j] * vet[lin[j]-1];
			}
			k++;
		}
	}		
	else{
		for(i=(id*(N/nproc));i<N;i++){ 
			res[k] = 0;
			for(j=ptr[i]-1;j<ptr[i+1]-1;j++){ 
				res[k] += val[j] * vet[lin[j]-1];
			}
			k++;
		}
	}
}

void produto_matriz_transposta_vetor(double *mat, int *lin, int *ptr, double *vet, double *res, int N){
	int i,j;
	// printf("---produto_matriz_transposta_vetor---\n");

	for(j=0;j<N;j++){
		res[j] = 0;
	}

	for(j=0;j<N;j++){
		for(i=ptr[j]-1;i<ptr[j+1]-1;i++){
			int index = lin[i]-1;
			res[index] += mat[i] * vet[j];
		}
	}		
}

void subtracao_vetor(int n, double *vetB, double *vaux, double *vetR){
	for(int i=0;i<n;i++){
		vetR[i] = vetB[i] - vaux[i];
	}
}

void copia_vetor(int n, double *vetOrigem, double *vetDestino){
	for(int i=0;i<n;i++){
		vetDestino[i] = vetOrigem[i];
	}
}

double produto_escalar_vetores(int n, double *vet1, double *vet2){
	int i;
	double prod=0;
	for(i=0;i<n;i++)
		prod += vet1[i] * vet2[i];
	return(prod);	
}

double escalar_vetor(int n, double escalar, double *vet, double *res){
	int i;
	for(i=0;i<n;i++)
		res[i] = escalar * vet[i];
}

void soma_vetor(int n, double *vet1, double *vet2, double *res){
	int i;
	for(i=0;i<n;i++)
		res[i] = vet1[i] + vet2[i];
}


int calculo_gatterv(int id,int nproc,int N,int *cont, int *deslocamento){
	int i, tam = 0,d = 0;
	
	for(i=0; i<nproc; i++){
		if ( i != nproc-1 ){
			cont[i] = N/nproc;
			deslocamento[i] = d;
			d += N/nproc;
		}
		else{
			cont[i] = N/nproc + N%nproc;
			deslocamento[i] = d;
			d += N/nproc + N%nproc;
		}
	}	
	
	if ( id != nproc-1 ){
		tam = N/nproc;
	}
	else{
		tam = N/nproc + N%nproc;
	}
	return tam;
}

void bigC(int id, int nproc, int naozeros, int n, double *matA, double *vetB, int *linhas, int *colptr, double *vetX, double *matA_t, int *colunas_t, int *linptr_t){
	
	double *vaux, *vaux_p, *vetR, *vetR2,*vetR2T, *vetP, *vetP2, rho, rho0, beta, *vetV, *vetV_p, alpha, *matT;
	int i,tam,*cont = NULL,*deslocamento= NULL;

	vaux = (double *) malloc (n*sizeof(double));
	vaux_p = (double *) malloc (n/nproc*sizeof(double));
	vetR = (double *) malloc (n*sizeof(double));
	vetR2 = (double *) malloc (n*sizeof(double));
	vetR2T= (double *) malloc (n*sizeof(double));
	vetP = (double *) malloc (n*sizeof(double));	
	vetP2 = (double *) malloc (n*sizeof(double));
	vetV = (double *) malloc (n*sizeof(double));
	vetV_p = (double *) malloc (n/nproc*sizeof(double));
	
	matT = (double *) malloc (n*n*sizeof(double));
	cont = (int *)malloc(nproc*sizeof(int));
	deslocamento = (int *)malloc(nproc*sizeof(int));
	// x = 0	
	inicializa_vetor(n, vetX, 0);

	i = 1;
	tam = calculo_gatterv(id,nproc,n,cont,deslocamento);
	
	// r = b - Ax
	inicializa_vetor(tam, vaux_p, 0);
	
	
	produto_matriz_vetor(id, nproc, matA_t,colunas_t,linptr_t,vetX,vaux_p,n);

	
	if(id ==0){
		escreve_vetorInt(id, cont, nproc);
	}	

	MPI_Allgatherv(vaux_p, tam, MPI_DOUBLE, vaux, cont, deslocamento, MPI_DOUBLE, MPI_COMM_WORLD); 
	if(id==0){
		printf("vaux\n");
		escreve_vetor(id, vaux,n);
	}
	
	
	inicializa_vetor(n, vetR, 0);
	subtracao_vetor(n, vetB, vaux, vetR);
	copia_vetor(n, vetR, vetR2);

	//p = p2 = 0
	inicializa_vetor(n, vetP, 0);
	inicializa_vetor(n, vetP2, 0);
	
	rho = 1;
	
	while(i < IMAX){
		rho0 = rho;
		
		//rho = r2 * r    
		rho = produto_escalar_vetores(n, vetR2,vetR);

		beta = rho / rho0;
		printf("id=%d;n=:%d;beta:%lf\n",id,n,beta);
		
		// p = r + beta * p
		escalar_vetor(n, beta,vetP, vaux);
		
		soma_vetor(n, vetR, vaux, vetP);
		
		// p2 = r2 + beta * p2
		escalar_vetor(n, beta, vetP2, vaux);
		soma_vetor(n, vetR2,vaux, vetP2);

		printf("vetP\n");
		escreve_vetor(id,vetP,n);
		// v = A*p        
		inicializa_vetor(tam, vetV_p, 0);
		produto_matriz_vetor(id, nproc, matA_t,colunas_t,linptr_t,vetP,vetV_p,n);
		MPI_Allgatherv(vetV_p, tam, MPI_DOUBLE, vetV, cont, deslocamento, MPI_DOUBLE, MPI_COMM_WORLD); 

		if(id==0){
			printf("vetV\n");
			escreve_vetor(id,vetV,n);
		}
		
		//alpha = rho / (p2*v)        	
		alpha = rho/produto_escalar_vetores(n, vetP2, vetV);
		//x = x + alpha * p
		escalar_vetor(n, alpha, vetP, vaux);
		soma_vetor(n, vetX, vaux, vetX);
		
	    // r * r
		if (produto_escalar_vetores(n, vetR,vetR) < ERRO * ERRO) 
		    break;
		    
		//r = r - alpha * v
        
		escalar_vetor(n,alpha, vetV, vaux);
	    subtracao_vetor(n, vetR, vaux, vetR);
		
        //r2 = r2 - alpha *A' * p2 
		inicializa_vetor(tam, vaux_p, 0);
		produto_matriz_vetor(id,nproc,matA,linhas,colptr,vetP2,vaux_p,n);
		MPI_Allgatherv(vaux_p, tam, MPI_DOUBLE, vaux, cont, deslocamento, MPI_DOUBLE, MPI_COMM_WORLD); 

		escalar_vetor(n, alpha, vaux, vaux);
		subtracao_vetor(n,vetR2, vaux, vetR2);
		i++;
	}
	printf("vetX:\n");
	escreve_vetor(id,vetX, n);
	printf("\n");
	printf("i:%d\n",i);
	/**/


	free(vaux);
	free(vaux_p);
	free(vetR);
	free(vetR2);
	free(vetR2T);
	free(vetP);
	free(vetP2);
	free(vetV);
	free(vetV_p);
	free(matT);
	free(cont);
	free(deslocamento);

}


/*-----------------------------------------------------------------*/

int main(int argc, char **argv){ //argv é um vetor de strings.  //argc = 2 porque tenho 2 parametros (nome e nome da matriz)

	double *mat = NULL, *mat_t = NULL;
	double *vetB, *vetX;
	double sigma_novo, sigma_0,alpha, sigma_velho,beta;
	double inicio, fim;
	int *linhas = NULL, *colptr = NULL, *colunas = NULL, *linptr = NULL;
	int M, N, naozeros,n,i,j,id, nproc;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	// printf("ID:%d\t NPROC:%d\n",id,nproc);

	if (argc != 2){
		printf ("%s < Arquivo HB >\n",argv[0]);
		exit(-1);
	}
	
	if(id == 0){
		le_matriz(argv[1], &M, &N, &naozeros, &colptr, &linhas, &mat);
	}

	MPI_Bcast(&naozeros,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&N,1,MPI_INT,0,MPI_COMM_WORLD);

	if(id!= 0){
		linhas  = (int *) malloc (naozeros * sizeof(int));
		colptr  = (int *) malloc ((N+1) * sizeof(int));
		mat = (double *) malloc (naozeros * sizeof(double));
	}


	MPI_Bcast(linhas,naozeros,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(colptr,N+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(mat,naozeros,MPI_DOUBLE,0,MPI_COMM_WORLD);

	escreve_matriz(id, M, N, naozeros,nproc, colptr, linhas, mat);

	//------------Calcular Transposta ------------
	mat_t = (double *) malloc (naozeros * sizeof(double));
	colunas  = (int *) malloc (naozeros * sizeof(int));
	linptr  = (int *) malloc ((N+1) * sizeof(int));

	if(id ==0){
		calcula_transposta(id,mat,linhas,colptr, mat_t, colunas, linptr, N, naozeros);
	}

	MPI_Bcast(colunas,naozeros,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(linptr,N+1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(mat_t,naozeros,MPI_DOUBLE,0,MPI_COMM_WORLD);

	//------------Inicializar vetB------------
	vetB = (double *) malloc (N*sizeof(double));
	
	if(id == 0)
		inicializa_vetor(N, vetB, 1);

	MPI_Bcast(vetB,N,MPI_DOUBLE,0,MPI_COMM_WORLD);

	
	vetX = (double *) malloc (N*sizeof(double));


	bigC(id, nproc, naozeros,N, mat, vetB,linhas, colptr, vetX, mat_t, colunas, linptr);
	/**/

	free(linhas);
	free(colptr);
	free(mat);
	free(mat_t);
	free(colunas);
	free(linptr);


	MPI_Finalize();
}
