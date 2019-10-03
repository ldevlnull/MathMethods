#include <stdio.h>
#include "math.h"

const static double LEFT = 0;
const static double RIGHT = 5;
const static int INIT_N = 10;
const static int SCALE = 20;

const static char* FILE_INPUT_PATH = "input.txt";
const static char* FILE_OUTPUT_PATH = "output.txt";
const static char* FILE_OUTPUT_WKX_PATH = "output_wkx.txt";
const static char* FILE_OUTPUT_ERROR_PATH = "output_error.txt";
const static char* FILE_FORMAT = "%le\t%le\n";

double f(double x){
    return sin(x);
}

int countNodes() {
    FILE *input = fopen(FILE_INPUT_PATH, "r");
    char ch;
    int n = 0;
    while ((ch = fgetc(input)) != EOF)
        if (ch == '\n')
            n++;
    fclose(input);
    return n;
}

void readNodes(double *x, double *y){
    FILE *input = fopen(FILE_INPUT_PATH, "rt");
    int i = 0;
    while(!feof(input)){
        fscanf(input, FILE_FORMAT, &x[i], &y[i]);
        i++;
    }
    fclose(input);
}

double wkx(int k, double X, double *x){
    double product = 1;
    for(int i = 0; i <= k; i++) {
        product *= (X-x[i]);
    }
    return product;
}
double dividedDifferences(int k, double *x, double *y){
    double sum = 0;
    for(int i = 0; i <= k; i++){
        double p = 1;
        for(int j = 0; j <= k; j++){
            if(j != i)
                p *= (x[i]-x[j]);
        }
        sum += y[i] / p;
    }
    return sum;
}
double interpolate(double X, int N, double *x, double *y){
    double sum = y[0];
    int k;
    for(k = 1; k <= N; k++){
        sum += wkx(k-1, X, x)*dividedDifferences(k, x, y);
    }
    return sum;
}

void tabulateInFile(){
    double step = (RIGHT-LEFT)/INIT_N;
    double x, y;
    FILE *input = fopen(FILE_INPUT_PATH, "wt");
    for(int i = 0; i <= INIT_N; i++){
        x = LEFT + i*step;
        y = f(x);
        fprintf(input, FILE_FORMAT, x, y);
    }
    fclose(input);
}

int main(){
    tabulateInFile();
    int N = countNodes() - 2;
    double x[1000], y[1000];
    readNodes(x, y);
    double average_error = 0,
            X = x[0],
            step = (x[N]-x[0]) / (SCALE*N);
    FILE *output, *output_wkx, *output_error;
    output = fopen(FILE_OUTPUT_PATH, "wt");
    output_wkx = fopen(FILE_OUTPUT_WKX_PATH, "wt");
    output_error = fopen(FILE_OUTPUT_ERROR_PATH, "wt");
    for(int i = 0; i <= SCALE*N; i++){
        double error = fabs(f(X) - interpolate(X, N, x, y));
        fprintf(output, "%lf\t%lf\t%lf\n", X, f(X), interpolate(X, N, x, y));
        fprintf(output_wkx, "%lf\t%.20lf\n", X, wkx(N, X, x));
        fprintf(output_error, "%lf\t%.15lf\n", X, error);
        X += step;
        average_error += error;
    }
    average_error /= N*SCALE;
    fprintf(output_error, "average error - %.20lf", average_error);
    fclose(output);
    fclose(output_wkx);
    fclose(output_error);
}
