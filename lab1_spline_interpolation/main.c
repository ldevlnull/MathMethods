//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>


//double func(double x) {
//    return 2*sin(x);
//}

//void saveInFile() {
//    FILE *fp = fopen(FILE_PATH, "w+");
//    for (double x = B_A; x <= B_B; x += 3) {
//        fprintf(fp, FORMAT, x, func(x), 3);
//    }
//    fclose(fp);
//}

//int countNodes() {
//    FILE *fp = fopen(FILE_PATH, "r");
//    char ch;
//    int N = 0;
//    while ((ch = fgetc(fp)) != EOF)
//        if (ch == '\n')
//            N++;
//    return N;
//}

//void readNodes(float *x, float *y, float *h, int N){
//    FILE *fp = fopen(FILE_PATH, "r");
//    for (int i = 0; i < N; i++) {
//        fscanf(fp, FORMAT, &x[i], &y[i]);
//    }
//    for(int i = 1; i < N; i++){
//        h[i] = x[i] - x[i-1];
//    }
//    fclose(fp);
//}

//void thomasAlgoritm(float *h, float *y, float *c, const int NODES_AMOUNT){
//    float A[NODES_AMOUNT], B[NODES_AMOUNT], C[NODES_AMOUNT], Y[NODES_AMOUNT], a[NODES_AMOUNT], b[NODES_AMOUNT];
//    A[0] = B[0] = C[0] = 0;
//    Y[0] = 1;
//    for(int i = 1; i < NODES_AMOUNT; i++){
//        A[i] = h[i-1];
//        B[i] = h[i];
//        C[i] = 2*(h[i-1]+h[i]);
//        Y[i] = 6*((y[i]-y[i-1])/h[i]-(y[i-1]-y[i-2])/h[i-1]);
//    }
//    A[NODES_AMOUNT-1] = 0;
//    a[0] = -A[0]/C[0];
//    b[0] = Y[0]/C[0];
//    for(int i = 1; i < NODES_AMOUNT-1; i++){
//        a[i] = - A[i]/(A[i]*a[i-1]+C[i]);
//        b[i] = (Y[i]-A[i]*b[i-1])/(A[i]*a[i-1]+C[i]);
//    }
//    c[NODES_AMOUNT-1] = (Y[NODES_AMOUNT-1] - A[NODES_AMOUNT-1]*b[NODES_AMOUNT-2])/(A[NODES_AMOUNT-1]*a[NODES_AMOUNT-2]+C[NODES_AMOUNT-1]);
//    for(int i = NODES_AMOUNT-1; i > 0; i--){
//        c[i-1] = a[i-1]*c[i]+b[i-1];
//    }
//    FILE *output = fopen("output.txt", "w+");
//    printf("test");
//    fprintf(output, "i\t a\t b\t c\t Y\t h\n");
//    for(int i = 0; i < NODES_AMOUNT; i++){
//       fprintf(output, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", i, a[i], b[i], c[i], Y[i], h[i]);
//       printf("%f\n", h[i]);
//    }
//    fclose(output);

//}

//void calculateCoefficients(float *a, float *b, float *c, float *d, float *x, float *y, float *h, const int NODES_AMOUNT) {
////    thomasAlgoritm(h, &y, &c, NODES_AMOUNT);

//    float A[NODES_AMOUNT], B[NODES_AMOUNT], C[NODES_AMOUNT], Y[NODES_AMOUNT], alfa[NODES_AMOUNT], beta[NODES_AMOUNT];
//    A[0] = B[0] = C[0] = 0;
//    Y[0] = 1;
//    for(int i = 1; i < NODES_AMOUNT; i++){
//        A[i] = h[i];
//        B[i] = h[i+1];
//        C[i] = 2*(h[i]+h[i+1]);
//        Y[i] = 6*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
//    }
//    A[NODES_AMOUNT-1] = 0;
//    alfa[0] = -A[1]/C[1];
//    beta[0] = Y[1]/C[1];
//    for(int i = 1; i < NODES_AMOUNT-1; i++){
//        alfa[i] = - A[i+1]/(A[i+1]*alfa[i-1]+C[i+1]);
//        beta[i] = (Y[i+1]-A[i+1]*beta[i-1])/(A[i+1]*alfa[i-1]+C[i+1]);
//    }
//    c[NODES_AMOUNT-1] = (Y[NODES_AMOUNT-1] - A[NODES_AMOUNT-1]*beta[NODES_AMOUNT-2])/(A[NODES_AMOUNT-1]*alfa[NODES_AMOUNT-2]+C[NODES_AMOUNT-1]);
//    for(int i = NODES_AMOUNT-1; i > 0; i--){
//        c[i-1] = alfa[i-1]*c[i]+beta[i-1];
//        printf("%f\t%f\n", alfa[i], beta[i]);
//    }
//    for(int i = 1; i < NODES_AMOUNT; i++){
//        a[i-1] = y[i];
//        d[i-1] = (c[i] - c[i-1])/h[i];
//        b[i-1] = (h[i]/2)*c[i] - (sqrt(h[i])/6)*d[i] + (y[i]-y[i-1])/h[i];
//    }
//    a[NODES_AMOUNT-1]=y[NODES_AMOUNT-1-1];
//    b[NODES_AMOUNT-1]=(y[NODES_AMOUNT-1]-y[NODES_AMOUNT-1-1])/h[NODES_AMOUNT-1]-(2.0/3.0)*h[NODES_AMOUNT-1]*c[NODES_AMOUNT-1];
//    d[NODES_AMOUNT-1]=-c[NODES_AMOUNT-1]/(3*h[NODES_AMOUNT-1]);
//}

//int main() {
////    saveInFile();
//    const int NODES_AMOUNT = countNodes();
//    float x[NODES_AMOUNT], y[NODES_AMOUNT], h[NODES_AMOUNT];
//    readNodes(&x, &y, &h, NODES_AMOUNT);
//    float a[NODES_AMOUNT], b[NODES_AMOUNT], c[NODES_AMOUNT], d[NODES_AMOUNT];
//    calculateCoefficients(&a, &b, &c, &d, &x, &y, &h, NODES_AMOUNT);
//    FILE *output = fopen("output.txt", "w+");
//    fprintf(output, "i\t\ta\t\tb\t\tc\t\td\n");
//    for(int i = 0; i < NODES_AMOUNT; i++){
//       fprintf(output, "%d\t%e\t%e\t%e\t%e\n", i, a[i], b[i], c[i], d[i]);
//    }
//    fclose(output);
//}
//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>

//const double leftBorder = 1;
//const double rightBorder = 10;
//const char *FILE_PATH = "input.txt";
//const char *FORMAT = "%f\t%f\t\n";
//const int STEP = 1;


//double f(double x) {
//    return 5*sin(x);
//}

//void saveInFile() {
//    FILE *fp = fopen(FILE_PATH, "w+");
//    for (double x = leftBorder; x <= rightBorder; x += STEP) {
//        fprintf(fp, FORMAT, x, f(x));
//    }
//    fclose(fp);
//}

//int countNodes() {
//    FILE *fp = fopen(FILE_PATH, "r");
//    char ch;
//    int N = 0;
//    while ((ch = fgetc(fp)) != EOF)
//        if (ch == '\n')
//            N++;
//    return N;
//}

//void readNodes(float *x, float *y, int N){
//    FILE *fp = fopen(FILE_PATH, "r");
//    for (int i = 0; i < N+1; i++) {
//        fscanf(fp, FORMAT, &x[i], &y[i]);
//    }
//    fclose(fp);
//}

//int main() {
////    saveInFile();
//    int n, i, j;
//    n = countNodes();
//    n--;
//    float x[n + 1], a[n + 1], h[n], A[n], l[n + 1], u[n + 1], z[n + 1], c[n + 1], b[n], d[n];
//    readNodes(&x, &a, n);

//    for (i = 0; i <= n - 1; ++i){
//        h[i] = x[i + 1] - x[i];
//    }

//    for (i = 1; i <= n - 1; ++i){
//        A[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];
//    }

//    l[0] = 1;
//    u[0] = 0;
//    z[0] = 0;

//    for (i = 1; i <= n - 1; ++i) {
//        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
//        u[i] = h[i] / l[i];
//        z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
//    }

//    l[n] = 1;
//    z[n] = 0;
//    c[n] = 0;

//    for (j = n - 1; j >= 0; --j) {
//        c[j] = z[j] - u[j] * c[j + 1];
//        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
//        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
//    }

//    printf("%2s %8s %8s %8s %8s\n", "i", "ai", "bi", "ci", "di");
//    for (i = 0; i < n; ++i){
//        printf("%2d %8.2f %8.2f %8.2f %8.2f\n", i, a[i], b[i], c[i], d[i]);
//    }
//    const int SCALE = 15;
//    float xx[n*SCALE], yy[n*SCALE];
//    for (i = 0; i <= SCALE*n; i++){
//        xx[i] = x[0]+i*0.1f;
//        yy[i] = f(xx[i]);
//    }
//    float s = 0, eps = 0;
//    FILE *output = fopen("output.txt", "w+");
//    fprintf(output, "p,P\t\tx\t\ty(x)\t\ts(x)\t\teps\n");
//    for (i = 0, j = 1; i <= SCALE*(n-1); i++){
//        s=a[j]+b[j]*(xx[i]-x[j-1])+c[j]*(xx[i]-x[j-1])*(xx[i]-x[j-1])+
//        d[j]*(xx[i]-x[j-1])*(xx[i]-x[j-1])*(xx[i]-x[j-1]);
//        eps = fabs(s-yy[i]);
//        fprintf(output,"%i,%i\t%f\t% f\t%f\t%f\n",i,j,xx[i],yy[i],s,eps);
//        if ((i!=0)&&(i%SCALE==0)) {
//            j++;
//        }
//    }
//    fclose(output);
//    return 0;
//}
#include<stdio.h>
#include<math.h>

void gaussEliminationLS(int m, int n, double a[m][n], double x[n-1]){
    int i,j,k;
    for(i=0;i<m-1;i++){
        /*//Partial Pivoting
        for(k=i+1;k<m;k++){
            //If diagonal element(absolute vallue) is smaller than any of the terms below it
            if(fabs(a[i][i])<fabs(a[k][i])){
                //Swap the rows
                for(j=0;j<n;j++){
                    double temp;
                    temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
            }
        }*/
        //Begin Gauss Elimination
        for(k=i+1;k<m;k++){
            double  term=a[k][i]/ a[i][i];
            for(j=0;j<n;j++){
                a[k][j]=a[k][j]-term*a[i][j];
            }
        }

    }
    //Begin Back-substitution
    for(i=m-1;i>=0;i--){
        x[i]=a[i][n-1];
        for(j=i+1;j<n-1;j++){
            x[i]=x[i]-a[i][j]*x[j];
        }
        x[i]=x[i]/a[i][i];
    }

}
void cSCoeffCalc(int n, double h[n], double sig[n+1], double y[n+1], double a[n], double b[n], double c[n], double d[n]){
    int i;
    for(i=0;i<n;i++){
        d[i]=y[i];
        b[i]=sig[i]/2.0;
        a[i]=(sig[i+1]-sig[i])/(h[i]*6.0);
        c[i]=(y[i+1]-y[i])/h[i]-h[i]*(2*sig[i]+sig[i+1])/6.0;
    }
}
void tridiagonalCubicSplineGen(int n, double h[n], double a[n-1][n], double y[n+1]){
    int i;
    for(i=0;i<n-1;i++){
        a[i][i]=2*(h[i]+h[i+1]);
    }
    for(i=0;i<n-2;i++){
        a[i][i+1]=h[i+1];
        a[i+1][i]=h[i+1];
    }
    for(i=1;i<n;i++){
        a[i-1][n-1]=(y[i+1]-y[i])*6/(double)h[i]-(y[i]-y[i-1])*6/(double)h[i-1];
    }
}
void printMatrix(int m, int n, double matrix[m][n]){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            printf("%lf\t",matrix[i][j]);
        }
        printf("\n");
    }
}
void copyMatrix(int m, int n, double matrix1[m][n], double matrix2[m][n]){
    int i,j;
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            matrix2[i][j]=matrix1[i][j];
        }
    }
}
main(){
    int m,i;
    printf("Enter the no. of data-points:\n");
    scanf("%d",&m);
    int n=m-1;  //Now (n+1) is the total no. of data-points, following our convention
    double x[n+1]; //array to store the x-axis points
    double y[n+1]; //array to store the y-axis points
    double h[n];   ////array to store the successive interval widths
    printf("Enter the x-axis values:\n");
    for(i=0;i<n+1;i++){
        scanf("%lf",&x[i]);
    }
    printf("Enter the y-axis values:\n");
    for(i=0;i<n+1;i++){
        scanf("%lf",&y[i]);
    }
    for(i=0;i<n;i++){
        h[i]=x[i+1]-x[i];
    }
    double a[n]; //array to store the ai's
    double b[n]; //array to store the bi's
    double c[n]; //array to store the ci's
    double d[n]; //array to store the di's
    double sig[n+1]; //array to store Si's
    double sigTemp[n-1]; //array to store the Si's except S0 and Sn
    sig[0]=0;
    sig[n]=0;
    double tri[n-1][n]; //matrix to store the tridiagonal system of equations that will solve for Si's
    tridiagonalCubicSplineGen(n,h,tri,y); //to initialize tri[n-1][n]
    printf("The tridiagonal system for the Natural spline is:\n\n");
    printMatrix(n-1,n,tri);
    //Perform Gauss Elimination
    gaussEliminationLS(n-1,n,tri,sigTemp);
    for(i=1;i<n;i++){
        sig[i]=sigTemp[i-1];
    }
    //Print the values of Si's
    for(i=0;i<n+1;i++){
        printf("\nSig[%d] = %lf\n",i,sig[i]);
    }
    //calculate the values of ai's, bi's, ci's, and di's
    cSCoeffCalc(n,h,sig,y,a,b,c,d);
    printf("The equations of cubic interpolation polynomials between the successive intervals are:\n\n");
    for(i=0;i<n;i++){
        printf("P%d(x) b/w [%lf,%lf] = %lf*(x-%lf)^3+%lf*(x-%lf)^2+%lf*(x-%lf)+%lf\n",i,x[i],x[i+1],a[i],x[i],b[i],x[i],c[i],x[i],d[i]);
    }



}
