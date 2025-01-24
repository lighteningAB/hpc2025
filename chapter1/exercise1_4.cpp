// recursive matrix product

#include <iostream>

//concatenate two square arrays of size (n,n)
double* horizontalcat(double* a, double* b, int n){
    double * returnable = new double [2*n*n];
    int index = 0;
    //changes which row
    for(int i = 0; i<n; i++){
        //switches from a to b
        for (int s = 0; s<2; s++){
        //iterates through a row
            for (int j = 0; j<n; j++){
                if(s==0){
                    *(returnable+index)=(a[i*n+j]);
                    index++;
                }
                else{
                    *(returnable+index)=(b[i*n+j]);
                    index++;
                }
            }
        }
    }
    return returnable;
}

//concatenate two rectangular arrays of size (2n, n)
double* vertcat(double* a, double* b, int n){
    double * returnable = new double [4*n*n];
    int index = 0;
    for(int i = 0; i<2*n*n; i++){
        *(returnable+index)=(a[i]);
        index++;
    }
    for(int i = 0; i<2*n*n; i++){
        *(returnable+index)=(b[i]);
        index++;
    }
    return returnable;
}

//add two arrays of size (n,n)
double* vecadd(double * a, double * b, int n){
    double * returnable = new double [n*n];
    for(int i = 0; i<n*n; i++){
        *(returnable+i)=*(a+i)+*(b+i);
    }
    return returnable;
}

double* strassen(double* a, double* b, int* positions, int subsize, int fullsize){
    // base case, outputs a 2x2 matrix, flattened
    if (subsize == 2) {
        // Implement 2x2 matrix multiplication here
        double* result = new double[4];
        double a0 = a[positions[0]];
        double a1 = a[positions[0]+1];
        double a2 = a[positions[0]+fullsize];
        double a3 = a[positions[0]+fullsize+1];
        double b0 = b[positions[0]];
        double b1 = b[positions[0]+1];
        double b2 = b[positions[0]+fullsize];
        double b3 = b[positions[0]+fullsize+1];
        result[0] = a0*b0 + a1*b2;
        result[1] = a0*b1 + a1*b3;
        result[2] = a2*b0 + a3*b2;
        result[3] = a2*b1 + a3*b3;
        return result;
    }
    // recursion
    else {
        // Calculate start positions for each quadrant
        int halfSize = subsize / 2;
        
        // Positions for first matrix A
        int a11_start = positions[0];
        int a12_start = positions[0] + halfSize;
        int a21_start = positions[0] + halfSize * fullsize;
        int a22_start = positions[0] + halfSize * fullsize + halfSize;
        
        // Positions for second matrix B (similar calculation)
        int b11_start = positions[1];
        int b12_start = positions[1] + halfSize;
        int b21_start = positions[1] + halfSize * fullsize;
        int b22_start = positions[1] + halfSize * fullsize + halfSize;
        
        // New positions array for recursive calls
        int ae[] = {a11_start, b11_start};
        int bg[] = {a12_start, b21_start};
        int af[] = {a11_start, b12_start};
        int bh[] = {a12_start, b22_start};
        int ce[] = {a21_start, b11_start};
        int dg[] = {a22_start, b21_start};
        int cf[] = {a21_start, b12_start};
        int dh[] = {a22_start, b22_start};
        
        double* ae_n = strassen(a, b, ae, halfSize, fullsize);
        double* bg_n = strassen(a, b, bg, halfSize, fullsize);
        double* af_n = strassen(a, b, af, halfSize, fullsize);
        double* bh_n = strassen(a, b, bh, halfSize, fullsize);
        double* ce_n = strassen(a, b, ce, halfSize, fullsize);
        double* dg_n = strassen(a, b, dg, halfSize, fullsize);
        double* cf_n = strassen(a, b, cf, halfSize, fullsize);
        double* dh_n = strassen(a, b, dh, halfSize, fullsize);
        // Recursive calls with updated positions
        return vertcat(
            horizontalcat(
                vecadd(ae_n, bg_n, halfSize), //a
                vecadd(af_n, bh_n, halfSize), //b
                halfSize
            ), 
            horizontalcat(
                vecadd(ce_n, dg_n, halfSize), //c
                vecadd(cf_n, dh_n, halfSize), //d
                halfSize
            ),
            halfSize
        );
    }
}


int main() {
    // Test matrices
    double matrix_a[64] = {
    2, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
};

double matrix_b[64] = {
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1
};

    // Initial positions for full matrices
    int initial_positions[] = {0, 0};

    // Perform Strassen multiplication
    double* result = strassen(matrix_a, matrix_b, initial_positions, 8, 8);

    // Print the result
    std::cout << "Resulting Matrix:\n";
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            std::cout << result[i * 8 + j] << " ";
        }
        std::cout << "\n";
    }

    // Free dynamically allocated memory
    delete[] result;

    return 0;
}
