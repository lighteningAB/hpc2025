// computer freq power spectrum from signal.txt, identify dominant frequencies, plot original and cleaned up version
#include <fstream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <iostream>
using namespace std;

int main()
{
    vector<double> signal;
    double temp;

    ifstream sig("signal.txt");
    while (!sig.eof())
    {
        sig >> temp;
        signal.push_back(temp);
    }
    sig.close();

    const int n = signal.size();
    const int n_out = n / 2 + 1;
    const double srate = 1000.0;

    fftw_complex *tmp = new fftw_complex[n_out];
    fftw_plan p = fftw_plan_dft_r2c_1d(n, signal.data(), tmp, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    double fbinsize = srate / n;

    std::vector<double> out(n_out, 0.0);
    for (int i = 0; i < n_out; ++i)
    {
        out[i] = sqrt((tmp[i][0] * tmp[i][0] + tmp[i][1] * tmp[i][1])) / (n);
        if (i > 0)
        {
            out[i] *= 2.0;
        }
        cout << i * fbinsize << "  " << out[i] << endl;
    }

    return 0;
}
