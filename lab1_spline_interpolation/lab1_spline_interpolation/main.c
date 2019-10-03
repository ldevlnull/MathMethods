#include <stdio.h>
#include <math.h>

static const int INIT_N = 100;
static const double LEFT = 0;
static const double RIGHT = 10;

double f(double x){
    return 2*sin(x/3);
}

void thomasAlgorithm(double *y, double *h, int N, double *c){
    double alfa[N], beta[N], hamma[N], delta[N], A[N], B[N];
    alfa[1] = hamma[1] = delta[1] = 0;
    beta[1] = 1;
    for(int i = 2; i < N; i++){
        alfa[i]=h[i-1];
        beta[i]=2*(h[i-1]+h[i]);
        hamma[i]=h[i];
        delta[i]=3*(((y[i]-y[i-1])/h[i])-((y[i-1]-y[i-2])/h[i-1]));
    }
    hamma[N-1] = 0;
    A[1] = -hamma[1] / beta[1];
    B[1] = delta[1] / beta[1];
    for(int i = 2; i < N-1; i++){
        A[i] = -hamma[i]/(alfa[i]*A[i-1]+beta[i]);
        B[i] = (delta[i]-alfa[i]*B[i-1])/(alfa[i]*A[i-1]+beta[i]);
    }
    c[N-1] = (delta[N-1]-alfa[N-1]*B[N-2])/(alfa[N-1]*A[N-2]+beta[N-1]);
    for(int i=N-1; i>1; i--){
        c[i-1]=A[i-1]*c[i]+B[i-1];
    }
}

void tabulateInFile(){
    int N = INIT_N;
    double x0 = LEFT, xn = RIGHT;
    double x[N], y[N], h[N];
    double step = (xn-x0)/(N);
    for (int i = 0; i <= N; i++){
        x[i] = x0 + (i)*step;
        h[i] = step;
        y[i] = f(x[i]);
    }
    FILE *finput = fopen("input.txt","w");
    for(int i = 0; i <= N;i++){
        fprintf(finput,"\t%le\t%le\t%le\n", x[i], y[i], h[i]);
    }
    fclose(finput);
}
int calculateNodes(){
    FILE *input = fopen("input.txt","r");
    char currentChar;
    int N = 0;
    while((currentChar = fgetc(input)) != EOF)
        if (currentChar == '\n')
            N++;
    fclose(input);
    return N;
}

void readDataFromFile(double *x, double *y, double *h, int N){
    FILE *input = fopen("input.txt","r");
    for(int i = 0; i < N; i++){
        fscanf(input,"%le\t%le\t%le\n", &x[i], &y[i], &h[i]);
    }
    fclose(input);
}

void calculateKoefficients(double *a, double *b, double *c, double *d, double *y, double *h, int N){
    for(int i = 1; i < N-1;i++){
        a[i] = y[i-1];
        b[i] = (y[i]-y[i-1])/h[i]-(h[i]/3)*(c[i+1]+2*c[i]);
        d[i] = (c[i+1]-c[i])/(3*h[i]);
    }
    a[N-1] = y[N-2];
    b[N-1] = (y[N-1]-y[N-2])/h[N-1]-(2.0/3.0)*h[N-1]*c[N-1];
    d[N-1] = -c[N-1]/(3*h[N-1]);
}

double error(double y, double s){
    return fabs(y - s);
}

double interpolate(double *a, double *b, double *c, double *d, double *x, double *xm, int i, int j){
    return a[j]+b[j]*(xm[i]-x[j-1])+c[j]*(xm[i]-x[j-1])*(xm[i]-x[j-1])+ d[j]*(xm[i]-x[j-1])*(xm[i]-x[j-1])*(xm[i]-x[j-1]);
}

int main(){
    tabulateInFile();
    int N = calculateNodes(), scale = 50;
    double x[N], y[N], h[N], xm[scale*N], ym[scale*N], a[N], b[N], c[N], d[N];
    readDataFromFile(x, y, h, N);
    thomasAlgorithm(y, h, N, c);
    calculateKoefficients(a, b, c, d, y, h, N);

    double step = h[0] / scale;
    for (int i = 0; i <= scale*(N-1); i++){
        xm[i] = x[0] + i*step;
        ym[i] = f(xm[i]);
    }
    FILE *output = fopen("output.txt","w");

    fprintf(output,"\tx\t\t\ty\t\t\ts\t\t\teps\n");

    double avgError = 0;
    for (int i = 0, j = 1; i <= scale*(N-2); i++){
        double res = interpolate(a, b, c, d, x, xm, i, j);
        fprintf(output,"%.15lf\t%.15lf\t%.15lf\t%lf\n", xm[i], ym[i], res, error(ym[i], res));
        avgError += error(ym[i], res);
        if (i!=0 && i%scale == 0)
            j++;
    }
    fprintf(output, "average error - %.15lf", avgError/(N*scale));
    fclose(output);
    return 0;
}
