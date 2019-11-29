#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const double LEFT = 0;
static const double RIGHT = 5;
static const int INIT_N = 30;
static const int SCALE = 20;

static const char* FILE_INPUT_PATH = "input.txt";
static const char* FILE_OUTPUT_PATH = "output.txt";
static const char* FILE_OUTPUT_A_PATH = "output_A.txt";
static const char* FILE_OUTPUT_B_PATH = "output_B.txt";
static const char* FILE_OUTPUT_C_PATH = "output_C.txt";
static const char* FILE_OUTPUT_ERROR_PATH = "output_error.txt";
static const char* FILE_FORMAT = "%lf\t%lf\n";

double f(double x) {
    return sin(x);
}

void tabulateInFile() {
    double h = (RIGHT-LEFT) / INIT_N;
    double x;
    FILE *input;
    input = fopen(FILE_INPUT_PATH, "wt");
    for(int i = 0; i <= INIT_N; i++) {
        x = LEFT + i*h;
        fprintf(input, FILE_FORMAT, x, f(x));
    }
    fclose(input);
}

int countNodes() {
    FILE *fdata = fopen(FILE_INPUT_PATH, "r");
    int N = 0;
    char one_char;
    while((one_char = fgetc(fdata)) != EOF) {
        if (one_char == '\n') N++;
    }
    return N-1;
}

void readNodes(double *X, double *Y, int n){
    FILE *input = fopen(FILE_INPUT_PATH, "r");
    for(int i = 0; i <= n; i++) {
        fscanf(input, FILE_FORMAT, &X[i], &Y[i]);
    }
    fclose(input);
}

void fillMatrixB(double **B, double *X, int n, int m) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j <= m; j++) {
            B[i][j] = 0;
            for(int k = 0; k <= n; k++)
                B[i][j] += pow(X[k], (i+j));
        }
    }
}

void fillMatrixC(double *C, double *X, double *F, int n, int m) {
    for(int i = 0; i < m; i++){
        C[i] = 0;
        for(int j = 0; j <= n; j++)
            C[i] += F[j] * pow(X[j], i);
    }
}

void gaussMethod(double **A, double *B, double *X, int m) {
    double max = 0, sum = 0;
    int maxIndex = 0;
    for(int i = 0; i < m-1; i++) {
        maxIndex = i;
        max = fabs(A[i][i]);
        for(int j = i + 1; j < m; j++) {
            if(fabs(A[j][i]) > max) {
                max = fabs(A[j][i]);
                maxIndex = j;
            }
        }
        if(maxIndex != i) {
            for(int j = i; j < m; j++) {
                A[i][j] += A[maxIndex][j];
                A[maxIndex][j] = A[i][j] - A[maxIndex][j];
                A[i][j] -= A[maxIndex][j];
            }
            B[i] += B[maxIndex];
            B[maxIndex] = B[i] - B[maxIndex];
            B[i] -= B[maxIndex];
        }
        for(int j = i+1; j < m; j++) {
            double t = A[j][i]/A[i][i];
            for(int k = i; k < m; k++)
                A[j][k] -= t*A[i][k];
            B[j] -= t*B[i];
        }
    }

    X[m-1] = B[m-1] / A[m-1][m-1];

    for(int i = m-1; i >= 0; i--) {
        sum = 0;
        for(int j = i + 1; j < m; j++) {
            sum += A[i][j]*X[j];
        }
        X[i] = (1/A[i][i])*(B[i] - sum);
    }
}

double approximate(double *A, double x, int m) {
    double sum = 0;
    for(int i = 0; i < m; i++) {
        sum += A[i] * pow(x, i);
    }
    return sum;
}

double dispersion(double *X, double *A, int n, int m) {
    double D = 0, sum = 0;
    for(int i = 0; i <= n; i++) {
        sum += pow((f(X[i]) - approximate(A, X[i], m-1)), 2);
    }
    D = sqrt(sum / (n+1));
    return D;
}

void initMatrix(double **B, double *C, double *X, double *Y, int m, int n){
    fillMatrixB(B, X, n, m);
    fillMatrixC(C, X, Y, n, m);

    FILE *outputB = fopen(FILE_OUTPUT_B_PATH, "w");
    for(int k = 0; k < m; k++) {
        for(int j = 0; j < m-1; j++)
            fprintf(outputB, "%lf\t", B[k][j]);
        fprintf(outputB, "%lf\n", B[k][m]);
    }
    fprintf(outputB, "================= END OF ITERATION, m = %i ========================", m);
    fclose(outputB);

    FILE *outputC = fopen(FILE_OUTPUT_C_PATH,"w");
    for(int j = 0; j < m; j++)
        fprintf(outputC, "%lf\n", C[j]);
    fprintf(outputC, "================= END OF ITERATION, m = %i ========================", m);
    fclose(outputC);
}

void printMatrixA(double *A, int m){
    FILE *outputA = fopen(FILE_OUTPUT_A_PATH, "w");
    for(int j = 0; j < m; j++)
        fprintf(outputA, "%lf\n", A[j]);
    fprintf(outputA, "================= END OF ITERATION, m = %i ========================", m);
    fclose(outputA);
}

int main() {
    tabulateInFile();

    int m = 0;
    printf("input m...\n");
    scanf("%i", &m);
        printf("%i\n", m);
        int n = countNodes();

        double *A = (double *)malloc(m * sizeof(double));
        double **B = (double **)malloc(m * sizeof(double *));
        double *C = (double *)malloc(m * sizeof(double));
        double *X = (double *)malloc((unsigned long)(n+1) * sizeof(double));
        double *Y = (double *)malloc((unsigned long)(n+1) * sizeof(double));

        for(int i = 0; i <= m; i++)
            B[i] = (double *)malloc(m * sizeof(double));

        readNodes(X, Y, n);
        initMatrix(B, C, X, Y, m, n);
        gaussMethod(B, C, A, m);
        printMatrixA(A, m);

        double h = (X[n]-X[0])/(SCALE*n), delta = 0, x = 0;
        FILE *outputError = fopen(FILE_OUTPUT_ERROR_PATH, "w");
        FILE *output = fopen(FILE_OUTPUT_PATH, "w");
        fprintf(outputError, "eps = %lf\n", dispersion(X, A, n, m));
        fprintf(output, "\tx\t\t\tf(x)\t\ta(x)\terror(x)\n");

        for(int i=0; i <= SCALE*n; i++) {
            x = X[0] + i*h;
            fprintf(output,"%lf\t%lf\t%lf\t%lf\n", x, f(x), approximate(A,x,m), fabs(f(x)-approximate(A,x,m)));
            delta += pow((f(x)-approximate(A,x,m)), 2);
        }
        fprintf(output, "================= END OF ITERATION, m = %i ========================", m);
        fclose(output);
        delta = sqrt(delta / (SCALE * n+1));
        fprintf(outputError, "delta = %lf\n", delta);
        fprintf(outputError, "================= END OF ITERATION, m = %i ========================", m);
        fclose(outputError);

    return 0;
}
