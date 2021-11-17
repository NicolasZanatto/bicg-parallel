	/*------------------------------------------------------*/
//BIBLIOTECAS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iohb.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

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
void escreve_matriz(int M, int N, int naozeros, int *colptr, int *linhas, double *valores){
	int i;
	
	printf("VALORES:\n");
	for(i=0; i<naozeros; i++){
		printf("%f ", valores[i]);
	}
	printf("\n\n");

	printf("LINHAS:\n");
	for(i=0; i<naozeros; i++){
		printf("%d ", linhas[i]);
	}
	printf("\n\n");

	printf("PTR:\n");
	for(i=0; i<M+1; i++){
		printf("%d ", colptr[i]);
	}
	printf("\n");

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

void escreveVetor(double *vet, int n){
	int i;
	for(i=0; i<n; i++){
		printf("% 5.2lf\n", vet[i]);
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
	//int mat[25] = {5, 0, 2, 0, 0, 2, 5, 0, 0, 0, 0, 2, 5, 0, 1, 1, 0, 0, 5, 2, 1, 0, 1, 2, 5};
	for( i=0; i<n*n; i++){
		r = rand() % 10;
		matA[i] = r;
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

void produto_matriz_vetor (double *mat, int *lin, int *ptr, double *vet, double *res, int N){
	int i,j;
	for(i=0;i<N;i++){
		res[i] = 0;
		for(j=ptr[i];j<ptr[i+1];j++){
			int index = lin[j];
			res[i] += mat[j] * vet[index];
		}
	}		
}

void produto_matriz_transposta_vetor(double *mat, int *lin, int *ptr, double *vet, double *res, int N){
	int i,j;
	for(j=0;j<N;j++){
		res[j] = 0;
		for(i=ptr[j];i<ptr[j+1];i++){
			int index = lin[i];
			res[index] += mat[i] * vet[j];
		}
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

void bigC(int naozeros, int n, double *matA, double *vetB, int *linhas, int *colptr){
	double *vetX, *vaux, *vetR, *vetR2,*vetR2T, *vetP, *vetP2, rho, rho0, beta, *vetV, alpha, *matT;
	int i;

	vetX = (double *) malloc (n*sizeof(double));
	vaux = (double *) malloc (n*sizeof(double));
	vetR = (double *) malloc (n*sizeof(double));
	vetR2 = (double *) malloc (n*sizeof(double));
	vetR2T= (double *) malloc (n*sizeof(double));
	vetP = (double *) malloc (n*sizeof(double));	
	vetP2 = (double *) malloc (n*sizeof(double));
	vetV = (double *) malloc (n*sizeof(double));
	
	matT = (double *) malloc (n*n*sizeof(double));

	
	// x = 0	
	inicializaVetor(n, vetX, 0);
	
    	i = 1;
	
	// r = b - Ax
	inicializaVetor(n, vaux, 0);
	produto_matriz_vetor(matA,linhas,colptr,vetX,vaux,n);
	// produtoMatrizVetor(n, matA, vetX, vaux);
	
	inicializaVetor(n, vetR, 0);
	subtracaoVetor(n, vetB, vaux, vetR);
	copiaVetor(n, vetR, vetR2);
	
	//p = p2 = 0
	inicializaVetor(n, vetP, 0);
	inicializaVetor(n, vetP2, 0);
	
	rho = 1;

	while(i < IMAX){
		rho0 = rho;
		
		//rho = r2 * r    
        	rho = produtoEscalarVetores(n, vetR2,vetR);

	        beta = rho / rho0;

                // p = r + beta * p
		escalarVetor(n, beta,vetP, vaux);
		somaVetor(n, vetR, vaux, vetP);
	        // p2 = r2 + beta * p2
		escalarVetor(n, beta, vetP2, vaux);
		somaVetor(n, vetR2,vaux, vetP2);
		// v = A*p        
		inicializaVetor(n, vetV, 0);
        	//produtoMatrizVetor(n, matA,vetP,vetV);
		produto_matriz_vetor(matA,linhas,colptr,vetP,vetV,n);
        	//alpha = rho / (p2*v)        	
		alpha = rho/produtoEscalarVetores(n, vetP2, vetV);

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
		//transpostaMatriz(n, matA, matT);

		produto_matriz_transposta_vetor(matA,linhas,colptr,vetP2,vaux,n);		
		//produtoMatrizVetor(n, matT,vetP2,vaux);	
	
		escalarVetor(n, alpha, vaux, vaux);
		subtracaoVetor(n,vetR2, vaux, vetR2);


		i++;
	}
	printf("vetX:\n");
	escreveVetor(vetX, n);
	printf("\n");
	printf("i:%d\n",i);
	
	
}


/*-----------------------------------------------------------------*/

int main(int argc, char **argv){ //argv é um vetor de strings.  //argc = 2 porque tenho 2 parametros (nome e nome da matriz)

	double *valores = NULL;
	double *vetB;
	double sigma_novo, sigma_0,alpha, sigma_velho,beta;
	double inicio, fim;
	int *linhas = NULL, *colptr = NULL;
	int M, N, naozeros,n,i,j,iteracao=0;
	vetB = (double *) malloc (n*sizeof(double));

	srand(time(NULL));   // Initialization, should only be called once.

	if (argc != 2){
		printf ("%s < Arquivo HB >\n",argv[0]);
		exit(-1);
	}

	le_matriz(argv[1], &M, &N, &naozeros, &colptr, &linhas, &valores);    	

	escreve_matriz(M, N, naozeros, colptr, linhas, valores);

	inicializaVetor(N, vetB, 1);
	escreveVetor(vetB,N);
	printf("NAOZEROS:%d\n",naozeros);

	bigC(naozeros,N, valores, vetB,linhas, colptr);
}
/*------------------------------------------------------*/
