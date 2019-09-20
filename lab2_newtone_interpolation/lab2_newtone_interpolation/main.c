#include <stdio.h>
#include <math.h>

static int n = 0;

const char* FILE_INPUT_PATH = "input.txt";
const char* FILE_OUTPUT_PATH = "output.txt";
const static char* FILE_FORMAT = "%lf\t%lf\n";

double f(double x){
    return 5*cos(x);
}

void saveInFile() {
    FILE *fp = fopen(FILE_INPUT_PATH, "w+");
    for (double x = 0; x <= 10; x += 2) {
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
    for(int i = 0; i < k; i++){
        double p = 1;
        for(int j = 0; j < k; j++)
            if(j!=1)
                p *= (x[i] - x[j]);
            sum += y[i] / p;
    }
    return sum;
}

double wkx(double X, double *x, double k){
    double product = 1;
    for(int i = 0; i < k; i++){
        product *= X - x[i];
    }
    return product;
}

double N(double X, double *x, double *y){
    double sum = y[0];
    for(int i = 1; i < n; i++){
        sum += wkx(X, y, i)*dividedDifferences(x, y, i);
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
    readNodes(&x, &y);
    FILE *output = fopen(FILE_OUTPUT_PATH, "w+");
    fprintf(output, "i\tx\tN(x)\ty(x)\teps\n");
    for(int i = 0; i < n; i++){
//        fprintf(output, "%i%lf\t%lf\t%lf\t%lf\t\n", i, x[i], N(x[i], &x, &y), y[i],defineError(x[i], &x, &y));
        fprintf(output, "%lf\t%lf\n", x[i], N(x[i], &x, &y));
        printf("%lf\t%lf\n", x[i], N(x[i], &x, &y));
    }
    fclose(output);
}
