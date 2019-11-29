#include <stdio.h>
#include <math.h>

static char* PATH_TABULATION = "tabulation.txt";
static char* PATH_OUTPUT = "output.txt";
static char* PATH_INPUT = "input.txt";
static double eps = 1e-10;
static const double MAX_ITERATIONS = 1e6;

double F(double x, double *a, double *b, int n){
    b[n] = a[n];
    for(int i = n-1; i >= 0; i--) {
        b[i] = a[i] + x * b[i+1];
    }
    return b[0];
}

double f(double x){
    return cos(x);
//    return pow(x, 3) - 3*pow(x, 2) - 13*x + 15;
//    return pow(x, 3) + 6*pow(x, 2) - 15*x + 22;
}

double derivative(double x){
    return -sin(x);
//    return 3*pow(x, 2) - 6*x - 13;
//    return 3*pow(x, 2) + 12*x - 15;
}

double derivative2(double x){
    return -cos(x);
//    return 6*x - 6;
//    return 6*x + 12;
}

void tabulateInFile(){
    double a = 0, b = 5, n = 50;
    double h = (b - a) / n;
    FILE *out = fopen(PATH_TABULATION, "w+");
    for(double x = a; x <= b; x += h){
        fprintf(out, "%lf\t%.20lf\n", x, (double)f(x));
    }
    fclose(out);
}

int readCoefficientsFromFile() {
    FILE* in = fopen(PATH_INPUT, "r");
    int N = 0;
    char one_char;
    while((one_char = fgetc(in)) != EOF) {
        if (one_char == '\n') ++N;
    }
    fclose(in);
    return N-1;
}

double fixedPointIteration(double x0) {
    double tau = -1 / derivative(x0);
    double x = x0;
    int iterations = 0;
    do {
        x0 = x;
        x = x0 + tau*f(x0);
        iterations++;
    } while(!(fabs(x - x0) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS);
    printf("Fixed Point Method Iterations - %d\n", iterations);
    if (iterations >= MAX_ITERATIONS)
        printf("\nFIXED POINT ITERATION REACHED ITERATIONS LIMIT!\n");
    return x;
}

double NewtonMethod(double x0) {
    double x = x0;
    int iterations = 0;
    do {
        x0 = x;
        x = x0 - f(x0) / derivative(x0);
        iterations++;
    } while (!(fabs(x - x0) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS);
    printf("Newton Method Iterations - %d\n", iterations);
    if (iterations >= MAX_ITERATIONS)
        printf("\nNEWTON METHOD REACHED ITERATIONS LIMIT!\n");
    return x;
}

double ChebyshevMethod(double x0) {
    double x = x0;
    int iterations = 0;
    do {
        x0 = x;
        x = x0 - f(x0) / derivative(x0) -0.5*pow(f(x0), 2) * derivative2(x0) / pow(derivative(x0), 3);
        iterations++;
    } while (!(fabs(x - x0) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS);
    if (iterations >= MAX_ITERATIONS)
        printf("\nCHEBYSHEV METHOD REACHED ITERATIONS LIMIT!\n");
    return x;
}

double SecantMethod(double x0) {
    double tau, x1;
    int iterations = 0;
    tau = - 1 / derivative(x0);
    x1 = x0 + tau * f(x0);
    double x = x1;
    x1 = x0;
    while (!(fabs(x - x1) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS) {
        x0 = x1;
        x1 = x;
        x = x1 - f(x1) * (x1-x0) / (f(x1) - f(x0));
        iterations++;
    }
    printf("Secant Method Iterations - %d\n", iterations);
    if (iterations >= MAX_ITERATIONS)
        printf("\nSECANT METHOD REACHED ITERATIONS LIMIT!\n");
    return x;
}

double dividedDifferences(double x0, double x1) {
    return (f(x1)-f(x0)) / (x1-x0);
}

double delta(double x0, double x1, double x2, int ind)
{
    double res = 0;
    double dd1 = dividedDifferences(x1, x2), dd2 = (dd1 - dividedDifferences(x0, x1)) / (x2 - x0),
            root = sqrt(pow((x2 - x1) * dd2 + dd1, 2) - 4 * dd2 * f(x2));
    if(ind == 0) {
        res = (-((x2-x1)*dd2+dd1) + root)/(2*dd2);
    } else if (ind == 1) {
        res = (-((x2-x1)*dd2+dd1) - root)/(2*dd2);
    }
    return res;
}

double ParabolaMethod(double x0) {
    double x1, x2, delta0, delta1, delta2, tau;
    int iterations = 0;
    tau = -1/ derivative(x0);
    x1 = x0 + tau*f(x0);
    x2 = x1 + tau*f(x1);
    double x = x2;
    x2 = x1;
    x1 = x0;
    while (!(fabs(x - x2) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS) {
        x0 = x1;
        x1 = x2;
        x2 = x;
        delta1 = delta(x0, x1, x2, 0);
        delta2 = delta(x0, x1, x2, 1);
        delta0 = fabs(delta1) < fabs(delta2) ? delta1 : delta2;
        x = x2+delta0;
        iterations++;
    }
    printf("Parabola Method Iterations - %d\n", iterations);
    if (iterations >= MAX_ITERATIONS)
        printf("\nPARABOLA METHOD REACHED ITERATIONS LIMIT!\n");
    return x;
}

double ReverseInterpolationMethod(double x0) {
    double x1, tau;
    int iterations = 0;
    tau = - 1 / derivative(x0);
    x1 = x0 + tau*f(x0);
    double x = x1;
    x1 = x0;
    while (!(fabs(x - x1) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS) {
        x0 = x1;
        x1 = x;
        x = -x0*f(x1)/(f(x0)-f(x1)) - x1*f(x0)/(f(x1)-f(x0));
        iterations++;
    }
    printf("Reverse Interpolation Method Iterations - %d\n", iterations);
    if (iterations >= MAX_ITERATIONS)
        printf("\nREVERSE INTERPOLATION METHOD REACHED ITERATIONS LIMIT!\n");
    return x;
}

double AitkenMethod(double x0) {
    double x1, x2 = 0, tau;
    tau = - 1 / derivative(x0);
    double x = x0;
    int iterations = 0;
    while (!(fabs(x - x2) <= eps) && !(fabs(f(x)) <= eps) && iterations < MAX_ITERATIONS) {
        x0 = x;
        x1 = x0 + tau*f(x0);
        x2 = x1 + tau*f(x1);
        x = x2 + pow(x2-x1, 2) / (2*x1-x2-x0);
        iterations++;
    }
    printf("Aitken Method Iterations - %d\n", iterations);
    if (iterations >= MAX_ITERATIONS)
        printf("\nAITKEN METHOD REACHED ITERATIONS LIMIT!\n");
    return x;
}

double NewtonGornerMethod(double *a, double *b, double *c, int n) {
    FILE* in = fopen(PATH_INPUT, "r");
    double x0 = 2.51, x = 0;

    for (int i = 0; i <= n; i++) {
        fscanf(in,"%le\n", &a[i]);
    }
    x = x0;
    int iterations = 0;
    do {
        x0 = x;
        b[n] = a[n];
        for(int i = n-1; i >= 0; i--) {
            b[i] = a[i] + x0 * b[i+1];
        }
        c[n] = b[n];
        for(int i = n-1; i >= 1; i--) {
            c[i] = b[i] + x0 * c[i+1];
        }
        x = x0 - b[0]/c[1];
        iterations++;
    } while (!(fabs(x-x0) <= eps) && !(fabs(F(x0, a, b, n)) <= eps) && iterations < MAX_ITERATIONS);
    printf("Newton-Gorner Method Point Iterations - %d\n", iterations);
    return x;
}

double LinaMethod(double *a, double *b, int n, FILE *output) {
    double alpha0 = 1, beta0 = 1, alpha1, beta1, p0, q0, p1, q1;
    alpha1 = alpha0;
    beta1 = beta0;
    int iterations = 0;
    do {
        alpha0 = alpha1;
        beta0 = beta1;
        p0 = -2 * alpha0;
        q0 = alpha0*alpha0 + beta0*beta0;
        b[n] = a[n];
        b[n-1] = a[n-1] - p0*b[n];
        for(int i = n-2; i >= 2; i--){
            b[i] = a[i] - p0*b[i+1] - q0*b[i+2];
        }
        q1 = a[0] / b[2];
        p1 = (a[1]*b[2] - a[0]*b[3])/(b[2]*b[2]);
        alpha1 = -p1 / 2;
        beta1 = sqrt(fabs(q1 - alpha1*alpha1));
        iterations++;
    }
    while (!(fabs(alpha1 - alpha0) <= eps) && !(fabs(beta1 - beta0) <= eps) && iterations < MAX_ITERATIONS);
    if(iterations >= MAX_ITERATIONS)
        printf("Lina Method REACHED ITERATIONS LIMIT\n");
    fprintf(output, "Result = %e + %ei\n", alpha1, beta1);
    printf("Method Lina Method Iterations - %d\n", iterations);
    return alpha1;
}

int main() {
    FILE *output = fopen(PATH_OUTPUT, "w");
    int n = readCoefficientsFromFile();
    double *a = (double *)malloc((n+1) * sizeof(double));
    double *b = (double *)malloc((n+1) * sizeof(double));
    double *c = (double *)malloc((n+1) * sizeof(double));
    double x0 = 5;
    fprintf(output, "\nFixed Point Iteration\n");
    fprintf(output, "Result = %e\n", fixedPointIteration(x0));
    fprintf(output, "\nNewton Method\n");
    fprintf(output, "Result = %e\n", NewtonMethod(x0));
    fprintf(output, "\nMethod_Chebysheva\n");
    fprintf(output, "Result = %e\n", ChebyshevMethod(x0));
    fprintf(output, "\nSecant Method\n");
    fprintf(output, "Result = %e\n", SecantMethod(x0));
    fprintf(output, "\nParabola Method\n");
    fprintf(output, "Result = %e\n", ParabolaMethod(x0));
    fprintf(output, "\nReverse Interpolation Method\n");
    fprintf(output, "Result = %e\n", ReverseInterpolationMethod(x0));
    fprintf(output, "\nAitken Method\n");
    fprintf(output, "Result = %e\n", AitkenMethod(x0));
    fprintf(output, "\nNewton-Gorner Scheme\n");
    fprintf(output, "Result = %e\n", NewtonGornerMethod(a, b, c, n));
    fprintf(output, "\nLina Method\n");
    LinaMethod(a, b, n, output);
    fclose(output);
    printf("FINISHED!\n");
}
