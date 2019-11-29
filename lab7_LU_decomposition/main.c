#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

static const int N = 100;
static unsigned long long size;

static const char* PATH_A = "matrixA.txt";
static const char* PATH_B = "matrixB.txt";
static const char* PATH_X = "matrixX.txt";
static const char* PATH_L = "matrixL.txt";
static const char* PATH_U = "matrixU.txt";

void generateMatrixesInFile(){
    FILE *matrix_A = fopen(PATH_A, "w");
    double a = 1, b = 9, r_numb = 0.0;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            fprintf(matrix_A,"%le\t", (a+b*(double)(rand()/1000)));
        }
        r_numb = a+b*(double)(rand() / RAND_MAX);
        fprintf(matrix_A,"%le\n",r_numb);
    }
    fclose(matrix_A);
    printf("Generating of Matrix A has been finished! Output file: %s\n", PATH_A);
}

void readMatrixFromFile(double **A){
    FILE *matrixA = fopen(PATH_A, "r+");
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j < N; j++) {
            fscanf(matrixA,"%le\t", &A[i][j]);
        }
        fscanf(matrixA,"%le\n",&A[i][N]);
    }
    fclose(matrixA);
}

void allocMemory(double **A, double **L, double **U){
    for (int i = 0; i <= N; i++) {
        A[i] = (double *)malloc((N+1) * sizeof(double));
        L[i] = (double *)malloc((N+1) * sizeof(double));
        U[i] = (double *)malloc(N * sizeof(double));
    }
}

void printMatrixes(double **L, double **U){
    FILE *matrixL = fopen(PATH_L, "w");
    FILE *matrixU = fopen(PATH_U, "w");
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j < N; j++) {
            fprintf(matrixL,"%.20lf\t",L[i][j]);
            fprintf(matrixU,"%.20lf\t",U[i][j]);
        }
        fprintf(matrixL,"%.20lf\n",L[i][N]);
        fprintf(matrixU,"%.20lf\n",U[i][N]);
    }
    fclose(matrixL);
    fclose(matrixU);
    printf("Matrix'es L and U were saved in files: %s and %s\n", PATH_L, PATH_U);
}

void printResultMatrix(double *X){
    FILE *matrixX = fopen(PATH_X, "w");
    for (int i = 1; i <= N; i ++) {
        fprintf(matrixX, "%f\n", X[i]);
    }
    fclose(matrixX);
    printf("Matrix X was saved in file: %s\n", PATH_X);
}

double error(double *x) {
    double max = 0;
    for (int i = 1; i <= N; i++) {
        if (fabs(x[i]) > max)
            max = fabs(x[i]);
    }
    return max;
}

double eps(double **A, double *B, double *X){
    double max = -99999;
    for(int i = 1; i <= N; i++){
        double sum = 0;
        for(int j = 1; j <= N; j++){
            sum += A[i][j] * X[j];
        }
        sum -= B[i];
        if(sum > max)
            max = sum;
    }
    return max;
}

void LU_decomposition(double **L, double **U, double **A) {
    double sum;
    for (int k = 1; k <= N; k++) {
        for (int i = k; i <= N; i++) {
            sum = 0;
            for (int j = 1; j <= k-1; j++) {
                sum+=L[i][j]*U[j][k];
            }
            L[i][k] = A[i][k] - sum;
        }
        for (int i= k+1 ; i<= N; i++) {
            sum = 0;
            for (int j = 1; j<= k-1; j++) {
                sum+=L[k][j]*U[j][i];
            }
            U[k][i]=(A[k][i]-sum)/L[k][k];
        }
    }
    printf("LU-decomposition completed!\n");
}

void LU_solve(double **L, double **U, double *B, double *X) {
    double Y[N+1], sum = 0;
    Y[1] = B[1] / L[1][1];
    for (int k=2; k <= N; k++) {
        sum = 0.0;
        for(int j = 1; j <= k-1; j++) {
            sum += L[k][j]*Y[j];
        }
        Y[k] = (B[k] - sum)/L[k][k];
    }
    X[N] = Y[N];
    for (int k = N-1; k >= 1; k--) {
        sum = 0.0;
        for(int j=k+1; j <= N; j++) {
            sum += U[k][j]*X[j];
        }
        X[k] = Y[k] - sum;
    }
}

void initB(double **A, double *B, double x0){
    FILE *matrixB = fopen(PATH_B, "w");
    for (int i = 1; i <= N; i++) {
        double sum = 0;
        for (int j = 1; j <= N; j++) {
            sum += A[i][j];
        }
        B[i] = sum * x0;
        fprintf(matrixB,"%le\n",B[i]);
    }
    fclose(matrixB);
    printf("Matrix B saved in file: %s\n", PATH_B);
}

void initLU(double **L, double **U){
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            L[i][j] = 0.0;
            U[i][j] = 0.0;
            if (i == j) {
                U[i][j] = 1.0;
            }
        }
    }
}

void generateB(double **A, double *X, double *B) {
    double sum;
    for(int i = 1; i <= N; i++) {
        sum = 0;
        for(int j=1; j <= N; j++) {
            sum += A[i][j]*X[j];
        }
        B[i] = sum;
    }
}

void accuracyAdaptation(double **A, double *X, double *B, double **L, double **U, double x0){
    int kmax = 1e+6;
    double r_eps = 0, eps = 1e-12;

    double *dX = (double *)malloc(size);
    double *R = (double *)malloc(size);
    double *B0 = (double *)malloc(size);

    for (int i = 1; i <= N; i ++) {
        dX[i] = X[i] - x0;
    }

    r_eps = error(dX);
    printf("Start = %.20f\n", r_eps);
    int k = 0;
    do {
        generateB(A, X, B0);
        for (int i = 1; i <= N; i++) {
            R[i] = B[i] - B0[i];
        }
        LU_solve(L, U, R, dX);
        for (int i = 1; i <= N; i++) {
            X[i] += dX[i];
        }
        if (k >= kmax) {
            printf("Reached limit of iteration!\n");
            break;
        }
        k++;
    } while((!(error(dX) < eps))&&(!(error(R) < eps)));
    memset(dX, 0, size);
    memset(B0, 0, size);
    memset(R, 0, size);
    printf("Iterations = %d\n", k);
}

int main() {
    size = ((N+1) * sizeof(double));
    generateMatrixesInFile();
    double x0 = 0.24;

    double *B = (double *)malloc((N+1) * sizeof(double));
    double *X = (double *)malloc((N+1) * sizeof(double));
    double **A = (double **)malloc((N+1) * sizeof(double *));
    double **L = (double **)malloc((N+1) * sizeof(double *));
    double **U = (double **)malloc((N+1) * sizeof(double *));

    allocMemory(A, L, U);
    readMatrixFromFile(A);
    initB(A, B, x0);
    initLU(L, U);
    LU_decomposition(L, U, A);
    printMatrixes(L, U);
    LU_solve(L, U, B, X);
    accuracyAdaptation(A, X, B, L, U, x0);

    printf("eps = %le\n", eps(A, B, X));
    memset(B, 0, size);
    memset(X, 0, size);
    memset(A, 0, size);
    memset(L, 0, size);
    memset(U, 0, size);
}
