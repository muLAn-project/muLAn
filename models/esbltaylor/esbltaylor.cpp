//==================================================================//
// This file is part of the project muLAn.
// Compute the point source binary lens amplification and the
// extended source binary-lens amplification using the taylor
// approximation.
//
// This code uses the following conventions:
// * lens consists in two bodies called L1 and L2;
// * the corresponding masses are called M1 and M2;
// * the distance between L1 and L2 is s > 0 in Einstein units;
// * M1 > M2;
// * the mass ratio is q = M2/M1;
// * the x-axis joins L1 and L2;
// * L1 has the coordinates (x=0, y=0);
// * L2 has the coordinates (x=-s, y=0).
//==================================================================//

#include <iostream>
#include <complex>
#include <vector>

#include "esbltaylor.h"

using namespace std;

//==================================================================//
// External routines declarations
//==================================================================//
extern"C"{
    void wrap_(double* roots_r, double* roots_i, int* first_3_roots_order_changed, double* poly_r, double* poly_i, int* polish_only);
}

//==================================================================//
// Functions
//==================================================================//
void polynome_binary_lens(double s, double q, complex<double> zeta, double *poly_r, double *poly_i){

    int i;
    complex<double> ii(0.0,1.0), *poly;

    poly = (complex<double> *) malloc(6 * sizeof(complex<double>));

    for(i=0; i<6; i++){
        poly_r[i] = poly[i].real();
        poly_i[i] = poly[i].imag();
    }

    poly[0] = -pow(s,2) * zeta;
    poly[1] = -s*(1+q)*((2+pow(s,2))*zeta + 2*s*pow(abs(zeta),2)) - pow(s,2) *q;
    poly[2] = -(1+q)*(s*q + pow(s,2) *(q-1)*conj(zeta) + (1+q+pow(s,2) *(2+q))*zeta + pow(abs(zeta),2) *(2*s*(2+q)+pow(s,2) *(1+q)*(s+conj(zeta))));
    poly[3] = (1+q)*(pow(s,2) * q - s*(1+q)*zeta + (2*s+pow(s,3) * (1+q) + pow(s,2) * (1+q)*conj(zeta))*conj(zeta) - 2*pow(abs(zeta),2) * (1+q)*(1+pow(s,2)+s*conj(zeta)));
    poly[4] = (1+q)*(s*(q-pow(abs(zeta),2) * (1+q)) + (1+q)*((1+2*pow(s,2)) - pow(abs(zeta),2) + 2*s*conj(zeta))* conj(zeta));
    poly[5] = pow(1+q,2) * (s+conj(zeta)) * conj(zeta);

    for(i=0; i<6; i++){
        poly_r[i] = poly[i].real();
        poly_i[i] = poly[i].imag();
    }

    free(poly);
}
//------------------------------------------------------------------//
complex<double> lens_equation(double s, double q, complex<double> z){
    complex<double> ii(0.0,1.0), zeta;
    zeta = z - 1.0 / (1.0 + q) * (1.0 / conj(z) + q / (conj(z) + s));
    return zeta;
}
//------------------------------------------------------------------//
complex<double> w_n(double s, double q, complex<double> z, int n){
    if (abs(z+s) < 1e-12) z = z + 1e-10;
    if (abs(z) < 1e-12) z = 1e-10;
    return (1.0 / (1.0+q)) * ( 1.0/pow(z,n) + q/pow(z+s,n) );
}
//------------------------------------------------------------------//
double amp_ps(double s, double q, complex<double> z){
    return 1.0/(1.0-pow(abs(w_n(s,q,z,2)),2));
}
//------------------------------------------------------------------//
double getamp_psbl(double s, double q, double* z_r, double* z_i, complex <double> zeta,
                   int* flag_im, int* whichreal){
    int j, n_real, whichparity[5];
    double amp_total, ldiff_lim, x, whatamp[5];
    complex<double> ii(0.0,1.0), zeta_test, images;

    // Initialization
    n_real = 0;
    amp_total = 0.0;
    ldiff_lim = 1e-3;
    *flag_im = 1;

    for (j=0; j<5; j++){
        // Find real images
        images = z_r[j] + ii*z_i[j];
        zeta_test = lens_equation(s, q, images);
        if(abs(zeta_test-zeta)<ldiff_lim){
            n_real++;
            whichreal[j] = 1;
            // Compute amplification
            x = amp_ps(s, q, images);
            whatamp[j] = abs(x);
            // Parity
            if(x>0) whichparity[j] = 1;
            else whichparity[j] = -1;
            // Increment total amplification
            amp_total += whatamp[j];
            }
        else{
            whichreal[j] = 0;
            whatamp[j] = 0.0;
            whichparity[j] = 0.0;
        }
    }
    if (n_real<3 | n_real>5){
        // cout << "Wrong number of images for y=" << zeta.imag() << endl;
        *flag_im = 0;
    }
    return amp_total;
}

//==================================================================//
vector< vector<double> > magnifcalc(vector< vector<double> > input_matrix){

    // ATTENTION Faire que rayshooting si pas bon nombre d'images, toujours, en 24.

    // Declarations
    // ------------
    //cout << "psbl.cpp" << endl;
    int i, j, n_dates, degree, first_3_roots_order_changed=0, polish_only=0;
    int flag_control_err, flag_control_im, whichreal[5];
    double pi=acos(-1.0), ldiff_lim, sigma, mu_center, a_half_plus, a_full_plus;
    double a_full_cross, a2rho2, a4rho4, mu_2, mu_4, s, q, rho, gamma, err;
    double derr[2], roots_r[5], roots_i[5], poly_r[6], poly_i[6];
    complex<double> ii(0.0,1.0), zeta, zeta0, zeta_test, images, voisins[12];
    vector <double> amplification, flag;
    vector< vector<double> > result;

    // Initializations
    // ---------------
    ldiff_lim = 0.0001;
    q = input_matrix[0][0];
    rho = input_matrix[0][1];
    gamma = input_matrix[0][2];
    degree = input_matrix[0][3];
    err = input_matrix[0][4];
    sigma = input_matrix[0][5];
    n_dates = input_matrix[1].size();

    // Loop over the dates
    // -------------------
    derr[0] = 0;
    derr[1] = 0;
    for (i=0; i<n_dates; i++){
        derr[0] = derr[1] - derr[0];
        // Source position
        zeta0 = input_matrix[2][i] + ii*input_matrix[3][i];
        s = input_matrix[1][i];

        // Point source
        // Coefficients of polynomial equation, degree 5
        polynome_binary_lens(s, q, zeta0, poly_r, poly_i);
        // Solve polynomial equation
        wrap_(roots_r, roots_i, &first_3_roots_order_changed, poly_r, poly_i, &polish_only);
        // Magnification PSBL
        flag_control_im = 1;
        flag_control_err = 1;
        mu_center = getamp_psbl(s, q, roots_r, roots_i, zeta0, &flag_control_im, whichreal);
        if(degree==0) amplification.push_back(mu_center);

        // Source size effects
        if(degree==2 || degree==4){
            voisins[0] = rho;
            voisins[1] = -rho*ii;
            voisins[2] = -rho;
            voisins[3] = rho*ii;
            a_full_plus = 0.0;
            a_half_plus = 0.0;
            // Quadrapole approximation
            for(j=0; j<4; j++){
                zeta = zeta0 + voisins[j];
                polynome_binary_lens(s, q, zeta, poly_r, poly_i);
                wrap_(roots_r, roots_i, &first_3_roots_order_changed, poly_r, poly_i, &polish_only);
                a_full_plus += getamp_psbl(s, q, roots_r, roots_i, zeta, &flag_control_im, whichreal);

                voisins[j+8] = 0.5*voisins[j];
                zeta = zeta0 + voisins[j+8];
                polynome_binary_lens(s, q, zeta, poly_r, poly_i);
                wrap_(roots_r, roots_i, &first_3_roots_order_changed, poly_r, poly_i, &polish_only);
                a_half_plus += getamp_psbl(s, q, roots_r, roots_i, zeta, &flag_control_im, whichreal);
            }
            a_full_plus = a_full_plus/4.0 - mu_center;
            a_half_plus = a_half_plus/4.0 - mu_center;
            a2rho2 = (16.0*a_half_plus - a_full_plus)/3.0;
            mu_2 = mu_center + 0.5*a2rho2*(1.0 - 0.2*gamma);
            derr[1] = abs(mu_2 - mu_center)/abs(mu_2) + derr[0];
            if(degree==2) amplification.push_back(mu_2);

            // Hexadecapole
            if(degree==4){
                a_full_cross = 0.0;
                for(j=0; j<4; j++){
                    voisins[j+4] = voisins[j]*exp(ii*pi/4.0);
                    zeta = zeta0 + voisins[j+4];
                    polynome_binary_lens(s, q, zeta, poly_r, poly_i);
                    wrap_(roots_r, roots_i, &first_3_roots_order_changed, poly_r, poly_i, &polish_only);
                    a_full_cross += getamp_psbl(s, q, roots_r, roots_i, zeta, &flag_control_im, whichreal);
                }
                a_full_cross = a_full_cross/4.0 - mu_center;
                a4rho4 = 0.5*(a_full_plus + a_full_cross) - a2rho2;
                mu_4 = mu_2 + a4rho4*(1.0 - gamma*11.0/35.0)/3.0;

                derr[1] = abs(mu_4 - mu_2)/abs(mu_4) + derr[0];
                amplification.push_back(mu_4);
                //cout << "Hex " << derr[1] << endl;
            }
        }
        if (derr[1] > err) flag_control_err = 0;
//        if (derr[1] > err) {
//            cout << degree << " " << derr[1] << " " << flag_control_err << endl;
//        }
        flag.push_back(flag_control_err * flag_control_im);
    }

    result.push_back(amplification);
    result.push_back(flag);

    return result;
}
