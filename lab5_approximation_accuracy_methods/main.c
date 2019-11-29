#include<stdio.h>
#include<math.h>

long double f(long double x){
    return sinl(x);
}

long double derivative(long double x){
    return cosl(x);
}

long double dividedDifferences(long double x, long double h){
    return (f(x+h)-f(x-h))/(2*h);
}

double eps(long double val1, long double val2){
    return (double)fabsl(val1 - val2);
}

long double ApproximationMethod(long double x, long double h){
    return dividedDifferences(x, h);
}

long double RungeRumbergMethod(long double x, long double h){
    return dividedDifferences(x, h)+(dividedDifferences(x, h)-dividedDifferences(x, 2.0L*h))/3.0L;
}

long double AitkenMethod(long double x, long double h){
    return (dividedDifferences(x, 2.0L*h)*dividedDifferences(x, 2.0L*h)-dividedDifferences(x, h)*dividedDifferences(x, 4L*h)) / (2.0L*dividedDifferences(x, 2.0L*h)-(dividedDifferences(x, h)+dividedDifferences(x, 4.0L*h)));
}

double EstimatedAccuracy(long double x, long double h){
    return (double)(1/logl(2)*logl((dividedDifferences(x, 4.0L*h)-dividedDifferences(x, 2.0L*h))/(dividedDifferences(x, 2.0L*h)-dividedDifferences(x, h))));
}

void printResult(FILE* fdata, long double x0, long double h0, long double h1){
    fprintf(fdata, "\n===== RESULTS =====\n\n");
    fprintf(fdata, "fp(x0) = %.20f\n", (double)derivative(x0));
    fprintf(fdata, "Optimal Step =%25.20f(%.0e)\n", (double)h0, (double)h0);
    fprintf(fdata, "Optimal Step =%25.20f (%.0e)\n", (double)h1, (double)h1);
    fprintf(fdata, "Error of Derivative Approximation with optimal step =%25.20f\n", eps(derivative(x0), ApproximationMethod(x0, h0)));
    fprintf(fdata, "Error of Derovatove Approximation with net step     =%25.20f\n", eps(derivative(x0), ApproximationMethod(x0, h1)));
    fprintf(fdata, "Error of Runge-Remberg Method Approximation         =%25.20f\n", eps(derivative(x0), RungeRumbergMethod(x0, h1)));
    fprintf(fdata, "Error of Aitken Method Approximation                =%25.20f\n", eps(derivative(x0), AitkenMethod(x0, h1)));
    fprintf(fdata, "Estimated Accuracy                                  =%25.20f\n", EstimatedAccuracy(x0, h1));
}

int main(){
    FILE *output = fopen("output.txt","w");
    long double x0 = 1l, h = 0.001l, h0 = 0l, h1;
    double min = 10.0l;

    fprintf(output, "====== START OF DERIVATIVE APPROXIMATION WITH H ======\n\n");
    for(int i = 20; i >= -3; i--){
        h = powl(10, -i);
        fprintf(output, "%25.20f\t%25.20f\n", (double)h, eps(derivative(x0), ApproximationMethod(x0, h)));
        if(eps(derivative(x0), ApproximationMethod(x0, h)) < min){
            min = eps(derivative(x0), ApproximationMethod(x0, h));
            h0 = h;
        }
    }
    fprintf(output, "\n====== END OF DERIVATIVE APPROXIMATION WITH H ======\n\n");
    printf("%f\n", (double)h0);
    long double por = log10l(h0), hh = 1;
    for(int i = 1; i <= fabsl(por-1); i++) {
        hh /= 10.0l;
    }

    fprintf(output, "===== START OF DERIVATIVE APPROXIMATION WITH 2H =====\n\n");
    for(long double h = hh; h <= 100.0l*hh; h += hh) {
        fprintf(output,"%25.20f\t%25.20f\n", (double)h, eps(derivative(x0), ApproximationMethod(x0, h)));
        if(eps(derivative(x0), ApproximationMethod(x0, h)) < min) {
            min = eps(derivative(x0), ApproximationMethod(x0, h));
            h0 = h;
        }
    }
    fprintf(output, "\n===== END OF DERIVATIVE APPROXIMATION WITH 2H =====\n\n");

    h1 = 1e-6l;
    printResult(output, x0, h0, h1);
    fclose(output);
    return 0;
}
