#include <math.h>
#include <stdio.h>

static const char* FILE_INPUT = "input.txt";
static const char* FILE_OUTPUT = "output.txt";
static const char* FILE_ERROR = "error.txt";
static const char* FILE_FORMAT = "%Lf\t%Lf\n";

static const long double LEFT = 0;
static const long double RIGHT = 1;
static const int SCALE = 20;
static const int INIT_N = 20;


long double f(long double x){
    return sinl(x);
}

void saveInFile(){
    long double step = 5;
    long double x = 0, y = 0;
    FILE *file1;
    file1 = fopen(FILE_INPUT, "wt");
    for(int i = 0; i <= INIT_N; i++) {
        x = LEFT + i*step;
        y = f(x);
        printf("%Lf", step);
        fprintf(file1, "%Le\t%Le\n", x, y);
    }
    fclose(file1);

}

void readNodes(long double *x, long double *y){
    FILE *input = fopen(FILE_INPUT, "rt");
    int i = 0;
    while(!feof(input)){
        fscanf(input, "%Le\t%Le\n", &x[i], &y[i]);
        i++;
    }
    fclose(input);
}

long double factorial(int n){
    if(n==0 || n==1)
        return 1;
    return n*factorial(n-1);
}

long double combination(int n, int k){
    return factorial(n) / (factorial(k) * factorial(n - k));
}

int step (int n){
    if (n%2)
    {
        return -1;
    }
    else
    {
        return 1;
    }

}

long double finiteDifference(int n, long double *y){
    long double r = 0;
    for(int k = 0; k <= n; k++){
        r += y[k]*step(n-k)*combination(n, k);
    }
    return r;
}

long double polynomial(long double t, int k) {
    long double p = 1;
    if(k == 0){
        p = 1;
    } else {
        for (int i = 0; i < k; i++){
            p *= (t-i);
        }
    }
    return p;
}

long double approximate(int n, long double t, long double *y) {
    long double res = 0;
    for(int k = 0; k <= n; k++){
        res += finiteDifference(k, y) * polynomial(t, k) / factorial(k);
    }
    return res;
}

long double error(long double f_value, long double approximation_value) {
    return fabsl(f_value - approximation_value);
}

int main(){
    saveInFile();
    int n = INIT_N;
    long double step = (LEFT-RIGHT) / INIT_N;
    long double x[n], y[n];
    readNodes(x, y);

    FILE *output = fopen(FILE_OUTPUT, "w"), *output_error = fopen(FILE_ERROR, "w");
    long double var = 0, var_step = (n-LEFT)*(SCALE*n);
    long double average_error = 0;
    for (int j = 0; j <= SCALE*n; j++) {
        long double err = error(f(LEFT + step*var), approximate(n, var, y));
        fprintf(output, "%Lf\t%Lf\n", var, approximate(n, var, y));
        fprintf(output_error, "%Lf\t%Lf\n", var, err);
        average_error += err;
        var += var_step;
    }
    average_error /= SCALE*INIT_N;
    printf("%.20Lf\n", average_error);
    fclose(output);
    fclose(output_error);
}
