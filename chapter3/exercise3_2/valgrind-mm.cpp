#include <cstdlib>
#include <ctime>
using namespace std;

int main() {
    int i;
    int j;
    int k;
    double f;
    constexpr int s = 1000;
    constexpr int S = s*s;

    double* A = new double[S];
    double* B = new double[s];
    double* C = new double[S];

    srand(time(0));

    for (i = 0; i < S; ++i) {
        A[i] = rand()*100.0/RAND_MAX;
        B[i] = rand()*100.0/RAND_MAX;
    }

    for (i = 0; i < s; ++i) {
        for (j = 0; j < s; ++j) {
            for (k = 0; k < s; ++k) {
                *(C++) += A[i*s+k] * B[k*s+j];
            }
        }
    }
    C -= s;

    delete[] A;
    delete[] B;

    return 0;
}
