// Copyright 2017 Clement Ranc
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.
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

#include "invraylocal.h"

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
double getamp_rayshooting(double xs, double ys, double s, double q,
    double gamma, double rho, double sigma, double size_rect_pix) {

    // Declarations
    int i, j, k, count, icount, nbinmax, flag;
    int nth, nl,il;
    int nxi,nyi,ngridi,im,ixi,iyi,ix,iy;
    int flag_im, whichreal[5];
    double pi=acos(-1);
    double amp, ampf, alim;
    double xmin, ymin;
    double ximin, ximax, yimin, yimax;
    double xip, yip;
    double th,thmin,thmax,dth;
    double x,y,xi,yi;
    double gridi,gridi2, sum,sumx,sumy,r2,r0,costh,br;
    double ximage[5], yimage[5];
    int *xibin, *yibin;
    double *mass, *xl, *yl;

    complex<double> ii(0.0,1.0), zeta, zeta0, zeta_test, images;
    double roots_r[5], roots_i[5], poly_r[6], poly_i[6];
    int first_3_roots_order_changed=0, polish_only=0;

    // Initialization
    nbinmax = 1000000;
    alim = 1.0;
    nth = 1000;
    nl=2;

    if (size_rect_pix<1) size_rect_pix = 1;

    xibin = (int *) malloc(nbinmax * sizeof(int));
    yibin = (int *) malloc(nbinmax * sizeof(int));

    mass = (double *) malloc(nl * sizeof(double));
    xl = (double *) malloc(nl * sizeof(double));
    yl = (double *) malloc(nl * sizeof(double));

    // Grid parameters
    gridi2 = 2.0*rho;

    ximin = -100.0;
    ximax =  100.0;
    yimin = -100.0;
    yimax =  100.0;
    nxi = int((ximax-ximin)/gridi2) + 1;
    ximax = ximin + gridi2*real(nxi);
    nyi = int((yimax-yimin)/gridi2) + 1;
    yimax = yimin + gridi2*real(nyi);

    thmin = 0.0;
    thmax = 2.0*pi;
    dth = (thmax-thmin)/real(nth);
    th = thmin - dth;

    x = xs + rho*cos(th);
    y = ys + rho*sin(th);
    zeta0 = x + ii * y;
    polynome_binary_lens(s, q, zeta0, poly_r, poly_i);
    wrap_(roots_r, roots_i, &first_3_roots_order_changed, poly_r, poly_i, &polish_only);
    amp = getamp_psbl(s, q, roots_r, roots_i, zeta0, &flag_im, whichreal);

    // Shooting grid
    icount = 0;

    while(th <= 2.0 * pi){
        dth = 2.0*pi/real(nth);
        if(amp >= 10.0) dth = 2.0*pi/(real(nth/10)*amp);
        if(amp >= 800.0) dth = 2.0*pi/40000.0;
        th = th + dth;
        x = xs + rho*cos(th);
        y = ys + rho*sin(th);
        zeta0 = x + ii * y;
        polynome_binary_lens(s, q, zeta0, poly_r, poly_i);
        wrap_(roots_r, roots_i, &first_3_roots_order_changed, poly_r, poly_i, &polish_only);
        amp = getamp_psbl(s, q, roots_r, roots_i, zeta0, &flag_im, whichreal);

        for(im=0; im<5; im++){
            if (whichreal[im]==1 || flag_im==0){
                ximage[im] = roots_r[im];
                yimage[im] = roots_i[im];
                xip = ximage[im];
                yip = yimage[im];
                ixi = int((xip-ximin)/gridi2)+1;
                iyi = int((yip-yimin)/gridi2)+1;
                for(i=-size_rect_pix; i<size_rect_pix+1; i++){
                    for(j=-size_rect_pix; j<size_rect_pix+1; j++){
                        flag=0;
                        for(k=0; k<icount; k++){
                            if(xibin[k]==ixi+i && yibin[k]==iyi+j) flag = 1;
                        }
                        if(flag==0){
                            icount = icount + 1;
                            if(icount > nbinmax){
                                cout << "nbinmax too small" << endl;
                                break;
                            }
                            xibin[icount] = ixi+i;
                            yibin[icount] = iyi+j;
                        }
                    }
                }
            }
        }
    }

    mass[0] = 1.0 / (1.0 + q);
    mass[1] = q / (1.0 + q);
    xl[0] = 0.0;
    yl[0] = 0.0;
    xl[1] = -s;
    yl[1] = 0.0;

    gridi = rho * pow((3.0*pi*pi/2.0) * (sigma*sigma) * alim, 1.0/3.0); // g: grid size.
    ngridi = int(gridi2/gridi) + 1;
    gridi = gridi2/real(ngridi);
    count = 0;
    sum = 0.0;

    //icount = 500;

    for(i=0; i<icount; i++){
        xmin = ximin + gridi2*real(xibin[i]-1);
        ymin = yimin + gridi2*real(yibin[i]-1);
        xi = xmin - gridi/2.0;
        for(ix=0; ix<ngridi; ix++){
            xi = xi + gridi;
            yi = ymin - gridi/2.0;

            for(iy=0; iy<ngridi; iy++){
                yi = yi + gridi;

                sumx = 0.0;
                sumy = 0.0;

                for(il=0; il<nl; il++){
                    r2=(xi-xl[il])*(xi-xl[il])+(yi-yl[il])*(yi-yl[il]);
                    sumx=sumx+(xi-xl[il])*mass[il]/r2;
                    sumy=sumy+(yi-yl[il])*mass[il]/r2;
                }
                x=xi-sumx;
                y=yi-sumy;
                r0 = sqrt(pow(x-xs,2)+pow(y-ys,2));
                if(r0 <= rho){
                    count = count + 1;
                    costh=sqrt(1.0 - pow(r0/rho,2));
                    br = 1.0 - gamma * (1.0 - 1.5 * costh);
                    sum = sum + br;
                }
            }

        }
    }

    ampf = sum * pow(gridi,2) / (pi*pow(rho,2));
    //cout << ampf << " " << sum << " " << gridi << " " << rho << endl;

    // ---

    // Free memory
    free(xibin);
    free(yibin);
    free(mass);
    free(xl);
    free(yl);

    // End of function
    return ampf;

}















//==================================================================//
vector< vector<double> > magnifcalc(vector< vector<double> > input_matrix){

    // Declarations
    // ------------
    int i, n_dates, rectangle_pixsize;
    double mu_r, s, q, rho, gamma, sigma;
    complex<double> ii(0.0,1.0), zeta, zeta0, zeta_test, images;
    vector <double> magnification, flag;
    vector< vector<double> > result;

    // Initializations
    // ---------------
    q = input_matrix[0][0];
    rho = input_matrix[0][1];
    gamma = input_matrix[0][2];
    sigma = input_matrix[0][5];
    rectangle_pixsize = input_matrix[0][6];
    n_dates = input_matrix[1].size();

    // Magnification for each date
    // ---------------------------
    for (i=0; i<n_dates; i++){
        // Source position
        zeta0 = input_matrix[2][i] + ii*input_matrix[3][i];
        s = input_matrix[1][i];
        // Magnification with finite source
        mu_r = getamp_rayshooting(zeta0.real(), zeta0.imag(), s, q, gamma, rho, sigma, rectangle_pixsize);
        magnification.push_back(mu_r);
        flag.push_back(1);
    }

    // Return the magnification and other information
    result.push_back(magnification);
    result.push_back(flag);

    return result;
}

















