#include <stdio.h>
#include <math.h>

static int n = 0;

const char* FILE_INPUT_PATH = "input.txt";
const char* FILE_OUTPUT_PATH = "output.txt";
const static char* FILE_FORMAT = "%lf\t%lf\n";

double f(double x){
    return 2*sin(x/3);
}

void saveInFile() {
    FILE *fp = fopen(FILE_INPUT_PATH, "w+");
    for (double x = 0; x <= 30; x+=1) {
        fprintf(fp, FILE_FORMAT, x, f(x));
    }
    fclose(fp);
}

void countNodes() {
    FILE *fp = fopen(FILE_INPUT_PATH, "r");
    char ch;
    while ((ch = fgetc(fp)) != EOF)
        if (ch == '\n')
            n++;
    fclose(fp);
}

void readNodes(double *x, double *y){
    FILE *fp = fopen(FILE_INPUT_PATH, "r");
    for (int i = 0; i < n; i++) {
        fscanf(fp, FILE_FORMAT, &x[i], &y[i]);
    }
    fclose(fp);
}

double dividedDifferences(double *x, double *y, double k){
    double sum = 0;
    for(int i = 0; i <= k; i++){
        double p = 1;
        for(int j = 0; j <= k; j++){
            if(j!=i){
                p *= (x[i] - x[j]);
            }
        }
        sum += y[i] / p;
    }
    return sum;
}

double wkx(double X, double *x, double k){
    double product = 1;
    for(int i = 0; i <= k; i++){
        product *= (X - x[i]);
    }
    return product;
}

double N(double X, double *x, double *y){
    double sum = y[0];
    for(int i = 1; i < n-1; i++){
        sum += wkx(X, x, i-1)*dividedDifferences(x, y, i);
    }
    return sum;
}

double defineError(double X, double *x, double *y){
    return fabs(f(X)-N(X, x, y));
}

int main(){
    saveInFile();
    countNodes();
    double x[n], y[n];
    readNodes(x, y);

    FILE *output = fopen(FILE_OUTPUT_PATH, "w+");
    fprintf(output, "i\t  x\t\t  N(x)\t\t  y(x)\t\t  eps\n");

    int tabulationAmount = 20;
    double var_x = x[0], step = (x[n-1]-x[0])/(tabulationAmount*n);
    for(double i = 0; i <= n*tabulationAmount; i++){
        fprintf(output, "%i\t%lf\t%lf\t%lf\t%lf\t\n", (int)(i), var_x, N(var_x, x, y), f(var_x), defineError(var_x, x, y));
        var_x += step;
    }
    fclose(output);
}
