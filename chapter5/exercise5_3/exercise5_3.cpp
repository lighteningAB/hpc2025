#include <iostream>
using namespace std;

#include <boost/gil/extension/io/jpeg.hpp>
using namespace boost::gil;


#define F77NAME(x) x##_

extern "C"
{
    void F77NAME(dgesvd)(const char &jobu, const char &jobvt, const int &m, const int &n,
                 double *a, const int &lda, double *s, double *u, const int &ldu, double *vt,
                 const int &ldvt, double* work, const int &lwork, int* info);
}

int main(){
    rgb8_image_t img;
    read_image("plane.jpg", img, jpeg_tag());

    const int n = img.width();
    const int m = img.height();

    cout << "Image dimensions: " << n << " x " << m << endl;

    // Extract the pixel data to three matrices
    double* r = new double[n*m];
    double* g = new double[n*m];
    double* b = new double[n*m];
    auto imgv = const_view(img);
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < m; y++) {
            r[x*m+y] = (int)at_c<0>(imgv(x,y));
            g[x*m+y] = (int)at_c<1>(imgv(x,y));
            b[x*m+y] = (int)at_c<2>(imgv(x,y));
        }
    }

    int numvec = min(m,n);
    double* s = new double[3*numvec];
    double* u = new double[3*m*m];
    double* v = new double[3*n*n];

    //calculate optimal workspace size
    double wkhelp;
    int info;
    F77NAME(dgesvd)('A','A', m, n, r, m, s, u, m, v, n, &wkhelp, -1, &info);

    //svd on three colors
    int lwork = (int)(wkhelp);
    double* wk = new double[lwork*3];
    F77NAME(dgesvd)('A','A', m, n, r, m, s, u, m, v, n, wk, lwork, &info);
    F77NAME(dgesvd)('A','A', m, n, g, m, s+numvec, u+m*m, m, v+n*n, n, wk+lwork, lwork, &info); //writing into empty space 
    F77NAME(dgesvd)('A','A', m, n, b, m, s+numvec*2, u+2*m*m, m, v+2*n*n, n, wk+lwork*2, lwork, &info); //writing into empty space 

    rgb8_image_t out(n, m);
    auto outv = view(out);

    int nvec = 250;

    //creating compressed image
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            double rr = 0;
            double gg = 0;
            double bb = 0;
            for (int k = 0; k < nvec; ++k) {
                rr += u[k*m+y]*v[x*n+k]*s[k];
                gg += u[m*m+k*m+y]*v[n*n+x*n+k]*s[numvec+k];
                bb += u[2*m*m+k*m+y]*v[2*n*n+x*n+k]*s[2*numvec+k];
            }
            outv(x,y) = rgb8_pixel_t(rr, gg, bb);
        }
    }

    // Write out compressed image
    write_view("output250.jpg", outv, jpeg_tag());

    nvec = 100;

    //creating compressed image
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            double rr = 0;
            double gg = 0;
            double bb = 0;
            for (int k = 0; k < nvec; ++k) {
                rr += u[k*m+y]*v[x*n+k]*s[k];
                gg += u[m*m+k*m+y]*v[n*n+x*n+k]*s[numvec+k];
                bb += u[2*m*m+k*m+y]*v[2*n*n+x*n+k]*s[2*numvec+k];
            }
            outv(x,y) = rgb8_pixel_t(rr, gg, bb);
        }
    }

    // Write out compressed image
    write_view("output100.jpg", outv, jpeg_tag());


    nvec = 50;

    //creating compressed image
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            double rr = 0;
            double gg = 0;
            double bb = 0;
            for (int k = 0; k < nvec; ++k) {
                rr += u[k*m+y]*v[x*n+k]*s[k];
                gg += u[m*m+k*m+y]*v[n*n+x*n+k]*s[numvec+k];
                bb += u[2*m*m+k*m+y]*v[2*n*n+x*n+k]*s[2*numvec+k];
            }
            outv(x,y) = rgb8_pixel_t(rr, gg, bb);
        }
    }

    // Write out compressed image
    write_view("output50.jpg", outv, jpeg_tag());

    nvec = 10;

    //creating compressed image
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            double rr = 0;
            double gg = 0;
            double bb = 0;
            for (int k = 0; k < nvec; ++k) {
                rr += u[k*m+y]*v[x*n+k]*s[k];
                gg += u[m*m+k*m+y]*v[n*n+x*n+k]*s[numvec+k];
                bb += u[2*m*m+k*m+y]*v[2*n*n+x*n+k]*s[2*numvec+k];
            }
            outv(x,y) = rgb8_pixel_t(rr, gg, bb);
        }
    }

    // Write out compressed image
    write_view("output10.jpg", outv, jpeg_tag());

    nvec = 5;

    //creating compressed image
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < m; ++y) {
            double rr = 0;
            double gg = 0;
            double bb = 0;
            for (int k = 0; k < nvec; ++k) {
                rr += u[k*m+y]*v[x*n+k]*s[k];
                gg += u[m*m+k*m+y]*v[n*n+x*n+k]*s[numvec+k];
                bb += u[2*m*m+k*m+y]*v[2*n*n+x*n+k]*s[2*numvec+k];
            }
            outv(x,y) = rgb8_pixel_t(rr, gg, bb);
        }
    }

    // Write out compressed image
    write_view("output5.jpg", outv, jpeg_tag());

    //deallocate memory
    delete[] r;
    delete[] g;
    delete[] b;
    delete[] s;
    delete[] u;
    delete[] v;
    delete[] wk;
}