#include <cmath>

void populateMatrix(double* M, const int N) {
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++j) {
            double alpha = pow(exp(i/N), 2.0);
            double beta  = pow(exp(j/N), 2.0);
            M[i*N+j] = alpha / beta;
        }
    }
}

void setup(const int N) {
    double* M = new double[N*N];
    populateMatrix(M, N);
    delete[] M;
}

int main() {
    const int N = 2000;

    setup(N);

    return 0;
}
