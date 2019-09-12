#include <stdio.h>
#include <stdlib.h>

const double A = -10;
const double B = 10;
const double STEP = 1;
const char *FILE_PATH = "output.txt";

double func(double x){
    return x*x-16;
}

void saveInFile(){
    FILE *fp = fopen(FILE_PATH, "w+");
    for (double x = A; x <= B; x+=STEP) {
        fprintf(fp, "f(%.2f)=%.2f\n", x, func(x));
    }
    fclose(fp);
}

void printFromFile(){
    FILE *fp = fopen(FILE_PATH, "r");
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char *result = malloc(fsize+1);
    fread(result, 1, fsize, fp);
    fclose(fp);
    printf(result);
}

int main(){
    saveInFile();
    printFromFile();
}
