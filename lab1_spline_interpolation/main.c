#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double A = 1;
const double B = 30;
const char *FILE_PATH = "input.txt";
const char *FORMAT = "%f %f\n";

double func(double x) {
    return 2*sin(x);
}

void saveInFile() {
    FILE *fp = fopen(FILE_PATH, "w+");
    for (double x = A; x <= B; x += 3) {
        fprintf(fp, FORMAT, x, func(x));
    }
    fclose(fp);
}

int countNodes() {
    FILE *fp = fopen(FILE_PATH, "r");
    char ch;
    int N = 0;
    while ((ch = fgetc(fp)) != EOF)
        if (ch == '\n')
            N++;
    return N;
}

void readNodes(float *x, float *y, int N){
    FILE *fp = fopen(FILE_PATH, "r");
    for (int i = 0; i < N; i++) {
        fscanf(fp, FORMAT, &x[i], &y[i]);
    }
    fclose(fp);
}

void calculateCoefficients(double *a, double *b, double *c, double *d, double *x, double *y, int NODES_AMOUNT) {
    for(int i = 0; i < NODES_AMOUNT; i++){
        a[i] = x[i];
    }
}

int main() {
//    saveInFile();
    const NODES_AMOUNT = countNodes();
    float x[NODES_AMOUNT], y[NODES_AMOUNT];
    readNodes(&x, &y, NODES_AMOUNT);
//    for(int i = 0; i < NODES_AMOUNT; i++){
//        printf("%f\t%f\n", x[i], y[i]);
//    }
    float a[NODES_AMOUNT], b[NODES_AMOUNT], c[NODES_AMOUNT], d[NODES_AMOUNT];
    calculateCoefficients(&a, &b, &c, &d, &x, &y, NODES_AMOUNT);
    for(int i = 0; i < NODES_AMOUNT; i++){
        printf("%f\n", a[i]);
    }
}
