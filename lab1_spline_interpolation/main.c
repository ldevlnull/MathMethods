#include <stdio.h>
#include <math.h>
#include <stdlib.h>

const double A = 0;
const double B = 15;
const char *FILE_PATH = "input.txt";
const char *FORMAT = "%e %e\n";

double func(double x) {
	return 2*sin(x);
}

void saveInFile() {
	FILE *fp = fopen(FILE_PATH, "w+");
	for (double x = A; x <= B; x += 1) {
		fprintf(fp, FORMAT, x, func(x));
	}
	fclose(fp);
}

void readFromFile(double *a, double *b, double *c, double *d, double *x, double *y) {
	FILE *fp = fopen(FILE_PATH, "r");
	char ch;
	int N = 0;
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == '\n')
			N++;
	}
	N--;
	 a = calloc(sizeof a, N+1);
	 x = calloc(sizeof x, N+1);
	double xx[20], yy[20];
	N--;
	rewind(fp);
	//for (int/* i = 0; i < NODES_AMOUNT; i++) {
	//	fscanf(fp, FORMAT, &xx[i], &yy[i]);
	//}*/
	fclose(fp);
	x = xx;
	y = yy;
}

void calculateCoefficients(double *a, double *b, double *c, double *d, double *x, double *y) {

}

int main() {
	saveInFile();
	/*double *a = NULL, *b = NULL, *c = NULL, *d = NULL, *x = NULL, *y = NULL;
	readFromFile(&a, &b, &c, &d, &x, &y);
	calculateCoefficients(&a, &b, &c, &d, &x, &y);*/

	int N = 0;
	FILE *data = fopen(FILE_PATH, "r");
	char ch;
	while ((ch = fgetc(data)) != EOF) {
		if (ch == '\n')
			N++;
	}
	double x[20], y[20];
	rewind(data);
	for (int i = 0; i < N; i++) {
		fscanf(data, FORMAT, &x[i], &y[i]);
	};
	for (int i = 0; i < N; i++) {
		printf("%e\n", x[i]);
	}

}
