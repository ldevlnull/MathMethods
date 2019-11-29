#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/*
    1 - Rosenbrock Function
    2 - Power Function
    3 - Root Function
    4 - Wood Function
    5 - Pauel Function
    6 - Miele Function
*/

static const int SELECTED_FUNCTION = 1;
static const char* OUTPUT_PATH = "output.txt";
static const int MAX_ITERATIONS = 1e6;
static const double X_CHANGE_INIT[] = { 1e-2, 1e-2, 1e-2, 1e-2 };

static const double X_INIT_DEFAULT[] = { -1,  2,  0,  0 };
static const double X_INIT_WOOD[]    = { -3, -1, -3, -1 };
static const double X_INIT_PAUEL[]   = { -3, -1,  0,  1 };
static const double X_INIT_MIELE[]   = {  1,  2,  2,  2 };

static double X_Change[] = { 1e-2, 1e-2, 1e-2, 1e-2 };

static int N = 4;
static double
       accuracy1 = 1e-10,
       accuracy2 = 1e-10,
       stepChange1 = 2,
       stepChange2 = 2;

double f(double *X) {
    switch(SELECTED_FUNCTION){
        case 1: return 100 * pow(X[0]*X[0] - X[1], 2) + pow(X[0] - 1, 2) ; // Rosenbrock Function
        case 2: return pow(10*pow((X[0] - X[1]), 2) + pow((X[0] -1), 2), 4); // Power Function
        case 3: return pow(10*pow((X[0] - X[1]), 2) + pow((X[0] -1), 2), 0.25); // Root Function
        case 4: return 100*(X[1] - X[0]*X[0])*(X[1] - X[0]*X[0]) +
                (1 - X[0])*(1 - X[0]) +
                90*(X[3] - X[2]*X[2])*(X[3] - X[2]*X[2]) +
                (1 - X[2])*(1 - X[2]) +
                10.1*((X[1] - 1)*(X[1] - 1) + (X[3] - 1)*(X[3] - 1)) +
                19.8*(X[1] - 1)*(X[3] - 1); // Wood Function
        case 5: return pow(X[0] + 10*X[1], 2) +
                5*pow(X[2] - X[3],    2) +
               10*pow(X[0] - X[3],    4) +
                  pow(X[1] - 2*X[2],  4); // Pauel Function
        case 6: return pow(exp(X[0] - X[1]), 4) + 100*pow((X[1] - X[2]), 6) +
                pow(tan(X[2] - X[3]), 4) + pow(X[0], 8) + pow((X[3] - 1), 2);
        default: printf("No function selected!"); return 0;
    }
}

bool equals(double *A, double *B) {
    double max = 0;
    for(int i = 0; i < N; i++)
        max = fabs(B[i] - A[i]) > max ? fabs (B[i] - A[i]) : max;
    return max < accuracy2;
}

void exploratoryMovesForPattern(double *X, double *X1) {
    for(int i = 0; i < N; i++)
        X1[i] = X[i];


    for (int i = 0; i< N; i++) {
        X1[i] = X[i] + X_Change[i];
        if (f(X1) < f(X))
            continue;
        else {
            X1[i] = X[i] - X_Change[i];
            if (f(X1) < f(X))
                continue;
            else
                X1[i] = X[i];
        }
    }
}

void patternSearch(double *X, double *X1, double *X2P) {
    for (int i = 0; i < N; i++)
        X2P[i] = X[i] + stepChange2 * (X1[i] - X[i]);
}

void exploratoryMoves(double *X, double *X1) {
    for(int i = 0; i < N; i++)
        X1[i] = X[i];

    for (int i = 0; i < N; i++)
        do {
            X1[i] = X[i] + X_Change[i];
            if (f(X1) < f(X))
                break;
            else {
                X1[i] = X[i] - X_Change[i];
                if (f(X1) < f(X))
                    break;
                else {
                    X_Change[i] /= stepChange1;
                    X1[i] = X[i];
                }
            }
        } while (X_Change[i] > accuracy1);
}

void initBasicPoint(double *X) {
    switch (SELECTED_FUNCTION) {
            case 1: case 2: case 3: N = 2; for(int i = 0; i < N; X[i] = X_INIT_DEFAULT[i], i++); return;
            case 4: N = 4; for(int i = 0; i < N; X[i] = X_INIT_WOOD[i], i++); return;
            case 5: N = 4; for(int i = 0; i < N; X[i] = X_INIT_PAUEL[i], i++); return;
            case 6: N = 4; for(int i = 0; i < N; X[i] = X_INIT_MIELE[i], i++); return;
    }
}

int main() {
    double
              X[] = { 0, 0, 0, 0 },
             X1[] = { 0, 0, 0, 0 },
             X2[] = { 0, 0, 0, 0 },
            X2P[] = { 0, 0, 0, 0 };
    initBasicPoint(X);

    FILE *output = fopen(OUTPUT_PATH, "w+");
    int iterations;
    for(iterations = 0; iterations < MAX_ITERATIONS; iterations++) {
        for (int i = 0; i < N; i++){
            X[i] = X1[i];
            X_Change[i] = X_CHANGE_INIT[i];
        }

        exploratoryMoves(X, X1);
        if (equals(X, X1))
            break;

        for (int i = 0; i < N; i++) {
            X2[i] = X1[i];
            X1[i] = X[i];
        }

        do {
            for (int i = 0; i < N; i++) {
                X[i] = X1[i];
                X1[i] = X2[i];
            }
            patternSearch(X, X1, X2P);
            for (int i = 0; i < N; i++)
                X_Change[i] = X_CHANGE_INIT[i];

            exploratoryMovesForPattern(X2P, X2);
        } while (!equals(X2P, X2) && f(X2) < f(X1));
    }

    fprintf(output, iterations >= MAX_ITERATIONS ? "Reached itearations limit!\n" : "Iterations = %d\n", iterations);
    fprintf(output, "X\n");

    for (int i = 0; i < N; i++)
        fprintf(output, "%.15f\n", X[i]);


    return 0;
}
