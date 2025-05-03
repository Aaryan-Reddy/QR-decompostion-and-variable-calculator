#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_matrix(int rows, int cols, double A[rows][cols]) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
            printf("%10.4f ", A[i][j]);
        printf("\n");
    }
    printf("\n");
}

void print_vector(int n, double v[n]) {
    for (int i = 0; i < n; i++)
        printf("%10.4f\n", v[i]);
    printf("\n");
}

void givens_qr(int m, int n, double A[m][n], double Q[m][m], double R[m][n]) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            Q[i][j] = (i == j) ? 1.0 : 0.0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            R[i][j] = A[i][j];
    for (int j = 0; j < n; j++) {
        for (int i = m - 1; i > j; i--) {
            double a = R[i - 1][j];
            double b = R[i][j];
            if (fabs(b) < 1e-10) continue;
            double r = sqrt(a * a + b * b);
            double c = a / r;
            double s = -b / r;
            for (int k = 0; k < n; k++) {
                double t1 = c * R[i - 1][k] - s * R[i][k];
                double t2 = s * R[i - 1][k] + c * R[i][k];
                R[i - 1][k] = t1;
                R[i][k] = t2;
            }
            for (int k = 0; k < m; k++) {
                double t1 = c * Q[k][i - 1] - s * Q[k][i];
                double t2 = s * Q[k][i - 1] + c * Q[k][i];
                Q[k][i - 1] = t1;
                Q[k][i] = t2;
            }
        }
    }
}

void multiply_matrix(int r1, int c1, int c2, double A[r1][c1], double B[c1][c2], double result[r1][c2]) {
    for (int i = 0; i < r1; i++)
        for (int j = 0; j < c2; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < c1; k++)
                result[i][j] += A[i][k] * B[k][j];
        }
}

void transpose(int rows, int cols, double A[rows][cols], double At[cols][rows]) {
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            At[j][i] = A[i][j];
}

int invert_matrix(int n, double A[n][n], double inverse[n][n]) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inverse[i][j] = (i == j) ? 1.0 : 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(A[i][i]) < 1e-10) return 0;
        double temp = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= temp;
            inverse[i][j] /= temp;
        }
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = 0; k < n; k++) {
                    A[j][k] -= A[i][k] * factor;
                    inverse[j][k] -= inverse[i][k] * factor;
                }
            }
        }
    }
    return 1;
}

int main() {
    int choice;
    printf("Choose an option:\n1. QR decomposition using Givens rotations\n2. Solve linear system using c = (AtA)^-1 At b\nEnter choice: ");
    scanf("%d", &choice);
    if (choice == 1) {
        int m, n;
        printf("Enter number of rows and columns of matrix A: ");
        scanf("%d %d", &m, &n);
        double A[m][n], Q[m][m], R[m][n];
        printf("Enter matrix A (row-wise):\n");
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                scanf("%lf", &A[i][j]);
        givens_qr(m, n, A, Q, R);
        printf("Matrix Q:\n"); print_matrix(m, m, Q);
        printf("Matrix R:\n"); print_matrix(m, n, R);
    } else if (choice == 2) {
        int eq, vars;
        printf("Enter number of equations and variables: ");
        scanf("%d %d", &eq, &vars);
        double A[eq][vars], At[vars][eq], AtA[vars][vars], invAtA[vars][vars];
        double b[eq], Atb[vars], c[vars];
        printf("Enter coefficients for each equation (total %d equations):\n", eq);
        for (int i = 0; i < eq; i++) {
            printf("Equation %d coefficients (%d vars): ", i+1, vars);
            for (int j = 0; j < vars; j++) scanf("%lf", &A[i][j]);
            printf("Enter RHS value: "); scanf("%lf", &b[i]);
        }
        transpose(eq, vars, A, At);
        multiply_matrix(vars, eq, vars, At, A, AtA);
        if (!invert_matrix(vars, AtA, invAtA)) {
            printf("Matrix is singular and cannot be inverted.\n");
            return 1;
        }
        for (int i = 0; i < vars; i++) {
            Atb[i] = 0.0;
            for (int j = 0; j < eq; j++) Atb[i] += At[i][j] * b[j];
        }
        for (int i = 0; i < vars; i++) {
            c[i] = 0.0;
            for (int j = 0; j < vars; j++) c[i] += invAtA[i][j] * Atb[j];
        }
        printf("Solution vector c:\n"); print_vector(vars, c);
    } else {
        printf("Invalid choice.\n");
    }
    return 0;
}

