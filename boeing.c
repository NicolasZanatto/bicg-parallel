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


void escreveMatriz(double *mat, int n){
    int i,j;
    for(i=0; i<n; i++){
	for(j=0; j<n; j++){
		printf("% 5.2lf\t", mat[i*n+j]);
    	}
    	printf("\n");
    }
}

void escreveVetor(int id, double *vet, int n){
	int i;
	printf("ID=%d\n",id);
	for(i=0; i<n; i++){
		printf("% 5.2lf\n", vet[i]);
	}
}

void escreveVetorInt(int id, int *vet, int n){
	int i;
	printf("ID=%d\n",id);
	for(i=0; i<n; i++){
		printf("%d\n", vet[i]);
	}
}

void inicializaVetor(int n, double *vet, int valor){
	for(int i=0; i<n; i++){
		vet[i] = valor;
	}
}

void inicializaMatriz(int n, double *matA){
	int i;
	double r;
	// int mat[25] = {5, 0, 2, 0, 0, 2, 5, 0, 0, 0, 0, 2, 5, 0, 1, 1, 0, 0, 5, 2, 1, 0, 1, 2, 5};
	int mat[25] = {4,1,0,0,5,0,4,1,0,0,0,0,4,1,0,0,0,0,4,1,0,0,0,0,4};
	for( i=0; i<n*n; i++){
		r = rand() % 10;
		matA[i] = mat[i];
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
	// val.size = 2 pq nproc = 5 
	// lin.size = 10
	// ptr.size = 6


	if(id != nproc-1){  //id=0;N=5;nproc=5
		for(i=(id*(N/nproc));i<(id+1)*(N/nproc);i++){
			// inic=0;fim=2 
			res[k] = 0;
			for(j=ptr[i]-1;j<ptr[i+1]-1;j++){

				// printf("ptr[i+id+1]:%d\n",ptr[i+id+1]);
				// printf("ID:%d;j:%d;val[%d]:%lf;lin[%d]=%d;vet[%d]:%lf\n",id,j,j,val[j],j,lin[j],j,vet[lin[j]]);
				res[k] += val[j] * vet[lin[j]-1];
				// printf("res[%d]:%lf\n\n",i,res[i]);
				
			}
			k++;
			// printf("i:%d;res:%lf\n",i,res[i]);
		}
	}		
	else{
		// printf("Entrou\n");
		// printf("id*(N/nproc):%d\n",(id*(N/nproc)));
		// printf("N=%d\n",N);
		for(i=(id*(N/nproc));i<N;i++){ 
			res[k] = 0;
			for(j=ptr[i]-1;j<ptr[i+1]-1;j++){ // id=0: 0 -- 3
													// id=1: 3 -- 5
													// id=2: 5 -- 7
													// id=3: 7 -- 9
													// id=4: 9 -- 10
				// printf("ID:%d;i=%d;j:%d;val[%d]:%lf;lin[%d]=%d;vet[%d]:%lf\n",id,i,j,j,val[j],j,lin[j],j,vet[lin[j]]);
				res[k] += val[j] * vet[lin[j]-1];
				// printf("res[%d]:%lf\n\n",i,res[i]);
				
			}
			k++;
			// printf("i:%d;res:%lf\n",i,res[i]);
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
			// printf("j=%d;i=%d;lin[%d]=%d;index=%d;mat[%d]=%lf;vet[%d]=%lf\n",j,i,i,lin[i],index,i,mat[i],j,vet[j]);
			res[index] += mat[i] * vet[j];
		}
		// printf("\n\n");
	}		
}

void subtracaoVetor(int n, double *vetB, double *vaux, double *vetR){
	for(int i=0;i<n;i++){
		vetR[i] = vetB[i] - vaux[i];
	}
}

void copiaVetor(int n, double *vetOrigem, double *vetDestino){
	for(int i=0;i<n;i++){
		vetDestino[i] = vetOrigem[i];
	}
}

double produtoEscalarVetores(int n, double *vet1, double *vet2){
	int i;
	double prod=0;
	for(i=0;i<n;i++)
		prod += vet1[i] * vet2[i];
	return(prod);	
}

double escalarVetor(int n, double escalar, double *vet, double *res){
	int i;
	for(i=0;i<n;i++)
		res[i] = escalar * vet[i];
}

void somaVetor(int n, double *vet1, double *vet2, double *res){
	int i;
	for(i=0;i<n;i++)
		res[i] = vet1[i] + vet2[i];
}


void transpostaMatriz(int n, double *mat, double *matT){
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
		    matT[j*n+i] = mat[i*n+j];
		}
}

void bigC(int id, int nproc, int naozeros, int n, double *matA, double *vetB, int *linhas, int *colptr, double *vetX){
	
	double *vaux, *vaux_p, *vetR, *vetR2,*vetR2T, *vetP, *vetP2, rho, rho0, beta, *vetV, *vetV_p, alpha, *matT;
	int i;

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

	
	// x = 0	
	inicializaVetor(n, vetX, 0);

	i = 1;
	
	// r = b - Ax
	inicializaVetor(n/nproc, vaux_p, 0);
	
	
	produto_matriz_vetor(id, nproc, matA,linhas,colptr,vetX,vaux_p,n);

	
	// printf("vaux_p\n");
	// escreveVetor(id, vaux_p,n/nproc);
	
	MPI_Allgather(vaux_p, n/nproc, MPI_DOUBLE, vaux, n/nproc, MPI_DOUBLE, MPI_COMM_WORLD);

	if(id==0){
		printf("vaux\n");
		escreveVetor(id, vaux,n);
	}
	
	
	// produtoMatrizVetor(n, matA, vetX, vaux);
	
	
	inicializaVetor(n, vetR, 0);
	subtracaoVetor(n, vetB, vaux, vetR);
	copiaVetor(n, vetR, vetR2);
	
	// printf("vetR2\n");
	// escreveVetor(vetR2,n);

	//p = p2 = 0
	inicializaVetor(n, vetP, 0);
	inicializaVetor(n, vetP2, 0);
	
	rho = 1;
	
	while(i < 2){
		rho0 = rho;
		
		//rho = r2 * r    
		rho = produtoEscalarVetores(n, vetR2,vetR);

		beta = rho / rho0;
		printf("id=%d;n=:%d;beta:%lf\n",id,n,beta);
		
                // p = r + beta * p
		escalarVetor(n, beta,vetP, vaux);
		
		somaVetor(n, vetR, vaux, vetP);
	    // printf("vetP\n");
		// escreveVetor(vetP,n);
		
		// p2 = r2 + beta * p2
		escalarVetor(n, beta, vetP2, vaux);
		somaVetor(n, vetR2,vaux, vetP2);

		printf("vetP\n");
		escreveVetor(id,vetP,n);
		//
		// v = A*p        
		inicializaVetor(n/nproc, vetV_p, 0);
		produto_matriz_vetor(id, nproc, matA,linhas,colptr,vetP,vetV_p,n);
		MPI_Allgather(vetV_p, n/nproc, MPI_DOUBLE, vetV, n/nproc, MPI_DOUBLE, MPI_COMM_WORLD);

		if(id==0){
			printf("vetV\n");
			escreveVetor(id,vetV,n);
		}
		/*
		//alpha = rho / (p2*v)        	
		alpha = rho/produtoEscalarVetores(n, vetP2, vetV);
		// printf("i=%d;alpha:%lf\n",i,alpha);
        	//x = x + alpha * p
		escalarVetor(n, alpha, vetP, vaux);
		somaVetor(n, vetX, vaux, vetX);
		
	    // r * r
		if (produtoEscalarVetores(n, vetR,vetR) < ERRO * ERRO) 
		    break;
		    
		//r = r - alpha * v
        
		escalarVetor(n,alpha, vetV, vaux);
	    subtracaoVetor(n, vetR, vaux, vetR);
		
        //r2 = r2 - alpha *A' * p2 
		produto_matriz_transposta_vetor(matA,linhas,colptr,vetP2,vaux,n);		
	
		// printf("produto-transposta-vaux\n");
		// escreveVetor(vaux,n);
		escalarVetor(n, alpha, vaux, vaux);
		subtracaoVetor(n,vetR2, vaux, vetR2);
		i++;
	}
	printf("vetX:\n");
	// escreveVetor(vetX, n);
	printf("\n");
	printf("i:%d\n",i);
	/**/
	i++;
	}
}


/*-----------------------------------------------------------------*/

int main(int argc, char **argv){ //argv é um vetor de strings.  //argc = 2 porque tenho 2 parametros (nome e nome da matriz)

	double *mat = NULL, *mat_p = NULL;
	double *vetB, *vetX;
	double sigma_novo, sigma_0,alpha, sigma_velho,beta;
	double inicio, fim;
	int *linhas = NULL, *colptr = NULL, *sendcounts, *displs;
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
	vetB = (double *) malloc (N*sizeof(double));
	
	if(id == 0)
		inicializaVetor(N, vetB, 1);

	MPI_Bcast(vetB,N,MPI_DOUBLE,0,MPI_COMM_WORLD);

	
	vetX = (double *) malloc (N*sizeof(double));


	bigC(id, nproc, naozeros,N, mat, vetB,linhas, colptr, vetX);
	/**/
	MPI_Finalize();
}
