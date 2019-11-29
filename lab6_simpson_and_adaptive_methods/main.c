#include <stdio.h>
#include <math.h>

static int iterationsCount = 0;

double f(double x) {
    return sin(x);
}

double integrate(double a, double b) {
    return -cos(b) + cos(a);
}

double SimpsonMethod(double a, double b, int N){
    double h = (b - a) / N;
    double res = f(a) + f(a + N*h);
    for(int i = 1; i < N; i++){
        res += i%2 == 0?2*f(a + i*h):4*f(a + i*h);
    }
    return res*(h/3);
}

double error(double ideal, double real){
    return fabs(ideal - real);
}

double AdaptiveAlgotithm(double a, double b) {
    double result = 0,
           h = (b - a) / 2,
           Integral1 = h/3*(f(a)+4*f(a+h)+f(b)),
           Integral2 = h/6*(f(a)+4*f(a+h/2)+f(a+h)) + h/6*(f(a+h)+4*f(a+(3/2)*h)+f(b)),
           delta = 1e-12;

    if (error(Integral2, Integral1) < delta) {
        result = Integral2;
    } else {
        iterationsCount++;
        result = AdaptiveAlgotithm(a, a+h) + AdaptiveAlgotithm(a+h, b);
    }
    return result;
}


int main() {
    int N = 24;
    double a = 0, b = 3;
    printf("math           [1; 3] =%26.20f\n", integrate(a, b));
    printf("Simpson Method [1; 3] =%26.20f\n", SimpsonMethod(a, b, N));
    printf("error with N=%d       =%26.20f\n", N, error(integrate(a, b), SimpsonMethod(a, b, N)));
    double min = error(integrate(a, b), SimpsonMethod(a, b, 10));
    int optimalN = 10;
    FILE *output = fopen("output.txt", "w");
    for(int i = 11; i < 1000; i++){
        fprintf(output, "%i\t%25.20f\n", i, error(integrate(a, b), SimpsonMethod(a, b, i)));
        if(min > error(integrate(a, b), SimpsonMethod(a, b, i))){
            min = error(integrate(a, b), SimpsonMethod(a, b, i));
            optimalN = i;
        }
    }
    printf("optimal N             =    %i\n", optimalN);
    printf("best accuracy with optiman N = %25.20f\n", error(integrate(a, b), SimpsonMethod(a, b, optimalN)));

    double RumbergIntegral = SimpsonMethod(a, b, N) + (SimpsonMethod(a, b, N) - SimpsonMethod(a, b, N/2)) / 15;
    double RumbergError = error(SimpsonMethod(a, b, N), RumbergIntegral);
    printf("rumberg       = %25.20f\n"
           "rubmerg error = %25.20f\n", RumbergIntegral, RumbergError);

    double AitkenIntegral = (pow(SimpsonMethod(a, b, N/2), 2) - SimpsonMethod(a, b, N)*SimpsonMethod(a, b, N/4))
            / (2*SimpsonMethod(a, b, N/2) - (SimpsonMethod(a, b, N) + SimpsonMethod(a, b, N/4)));
    double AitkenError = error(SimpsonMethod(a, b, N), AitkenIntegral);
    printf("aitken        = %25.20f\n"
           "aitken error  = %25.20f\n", AitkenIntegral, AitkenError);

    printf("adaptive algotihm        = %25.20f\n", AdaptiveAlgotithm(a, b));
    printf("adaptive algorithm error = %25.20f\n", error(SimpsonMethod(a, b, optimalN), AdaptiveAlgotithm(a, b)));
    printf("iterations = %i\n\n", iterationsCount);
    return 0;
}
