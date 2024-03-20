/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   main.cpp
  * Author: Никита
  *
  * Created on 30 сентября 2021 г., 22:36
  */

#include <cstdlib>
#include<iostream>
#include<math.h>
#include <complex>
#include<fstream>


using namespace std;

/*
 *
 */
const int n = 20;
const int n2 = n * n;
const double pi = acos(-1);
const complex<double> icomp(0, 1);
double Tera = pow(10, 12);
double Giga = pow(10, 9);
double Mega = pow(10, 3);
double Hz = 3.0 * Giga;
double ewave = 299792456; // m/sec;
const complex<double> k0((2.0 * pi * Hz) / ewave, 0), k1 = 1.5 * k0;

const double lymda = abs(2.0 * pi / k0.real());
//const double lymda =1;
const double a = -lymda/2;
const double b = lymda / 2;
const double c = -lymda / 2;
const double d = lymda / 2;



const double h1 = (b - a) / n;
const double h2 = (d - c) / n;

complex<double>  Green(double p) {
    return ( 0.25) *icomp* (_j0(p)* _y0(p));
    //return(1.0 / (4.0 * icomp) * exp(icomp * p));
    //return(exp(icomp * k0 * p) / (4 * pi * p));

}


complex<double>  Ker(double x1, double y1, double x2, double y2) {
    complex<double> zero(0, 0);
    double p = sqrt(pow(x1 - x2 , 2) + pow(y1 - y2 , 2));
    if (p > 1e-10) {
        return Green(p);
    }
    else return(zero);
}
complex<double>  K(double x1, double y1) {
    //double p = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
    return k1 * k1*(1.0+0*abs(x1-y1)) - k0 * k0;
}


//double middlepryam(double a, double b, double xi) {
//    double nn = 1000, h, x, in = 0, i;
//    h = (b - a) / nn;
//    x = a + (h / 2);
//    while (x <= b - (h / 2)) {
//        in = in + (K(xi, x) * h);
//        x = x + h;
//    }
//    return in;
//}

//complex<double> middlepryam1(double x1, double y1, double a1, double b1, complex<double> phi[n], int ires, int num) {
//    double nn = 100, h1, t, t2, i, c, h2;
//    complex<double>in(0, 0);
//    h1 = (b1 - a1) / nn;
//    for (i = 0; i < nn; i++) {
//        t = a1 + (i + 0.5) * h1;
//        in += (K(x1, y1, x(t, num), xii(t, num)) * phi[ires] * (sqrt(prx(t, num) * prx(t, num) + pry(t, num) * pry(t, num)))) * h1;
//    }
//    return in;
//}
complex<double> middlepryam2(double a1, double b1, double a2, double b2, double xi1, double xi2) {
    double nn = 8, h22, h11, x2, x1, i = 0, c, t1, t2, y1,y2;
    complex<double> in(0, 0);

    h22 = (b2 - a2) / nn;
    h11 = (b1 - a1) / nn;

    for (int i1 = 0; i1 < nn; i1++) {
        y1 = a1 + (i1 + 0.5) * h11;
        for (int i2 = 0; i2 < nn; i2++) {
            y2 = a2 + (i2 + 0.5) * h22;
             in += Ker(y1,y2,xi1,xi2) * h11 * h22*K(y1,y2);
        }
    }
    return in;
}
double del(int i, int j) {

    if (i == j) {
        return(1);
    }
    else
        return(0);


}
//double del2(int i1, int j1, int i2, int j2) {
//    double  res = (i1 == i2 && j1 == j2);
//     
//    return res;
//}


double phi(double x1,double x2, int i, int j, double massx1[n], double massx2[n]) {
    if ((x1 >= massx1[i]) && (x1 <= massx1[i + 1]) && (x2 >= massx2[j]) && (x2 <= massx2[j + 1])) {
        return 1;
    }
    else {
        return 0;
    }

}




complex<double> u(double xi1, double xi2, complex<double> cc[n],double massx1[n],double massx2[n]) {
    int i,j;
    complex<double>s = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
             s = s + cc[i*n+j] * phi(xi1,xi2 ,i,j, massx1, massx2);
        }
       

    }
    return(s);



}

void Gauss(int k, complex<double> Matrix[n*n][n*n+ 1]) {
    
    if (Matrix[k][k] != (1.0, 0.0)) {
        complex<double> T = Matrix[k][k];
        for (int j = k; j < n2 + 1; j++) {
            Matrix[k][j] = Matrix[k][j] / T;
        }
    }
    for (int i = 0; i < n*n; i++) {
        if ((Matrix[i][k] != complex<double>(0.0, 0.0)) && (i != k)) {
            complex<double> T = Matrix[i][k];
            Matrix[i][k] = complex<double>(0.0, 0.0);
            for (int j = k + 1; j <n*n+1; j++) {
                Matrix[i][j] -= Matrix[k][j] * T;
            }
        }
    }
    if (k < n2 - 1) {
        Gauss(k + 1, Matrix);
    }
}

void GaussInv(complex<double> A[n * n][n * n + 1], complex<double>  bigA[n*n][2*n*n]) {
    complex<double>ed[n][n], bigA[n][2 * n], K, res[n][n];
    int i, j, k;
    for (i = 0; i < n; i++) {
        ed[i][i] = 1;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            bigA[i][j] = A[i][j];
            bigA[i][j + n] = ed[i][j];
        }
    }
    for (k = 0; k < n; k++) {
        for (i = 0; i < 2 * n; i++) {
            bigA[k][i] = bigA[k][i] / A[k][k];
        }
        for (i = k + 1; i < n; i++) {
            K = bigA[i][k] / bigA[k][k];
            for (j = 0; j < 2 * n; j++) {
                bigA[i][j] = bigA[i][j] - bigA[k][j] * K;
            }
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A[i][j] = bigA[i][j];
            }
        }
    }
    for (k = n - 1; k > -1; k--) {
        for (i = 2 * n - 1; i > -1; i--) {
            bigA[k][i] = bigA[k][i] / A[k][k];
        }
        for (i = k - 1; i > -1; i--) {
            K = bigA[i][k] / bigA[k][k];
            for (j = 2 * n - 1; j > -1; j++) {
                bigA[i][j] = bigA[i][j] - bigA[k][j] * K;
            }
        }
    }

    //for(i=0;i<n;i++){
    //    for (j = 0; j < n; j++) {
    //        res[i][j] = bigA[i][j + n];
    //    }
    //}
}

int main(int argc, char** argv) {
    double h ,x1[n + 1], x2[n + 1];//xi[n * n][2] 
    complex<double> A[n * n][n * n + 1], cc[n * n];
    int i, j, k;
    ofstream  out1("1ecr.txt");

   // h = (b - a) / n;



    for (i = 0; i < n + 1; i++) {

        x1[i] = a + i * h1;
        cout << x1[i] << endl;
    }
    for (i = 0; i < n + 1; i++) {

        x2[i] = c + i * h2;
        cout << x2[i] << endl;
    }
    //for (i = 0; i < n; i++) {
    //    for (j = 0; j < n; j++) {
    //         xi[i*n+j][0] = x1[j] + (h1 / 2);
    //         xi[i* n + j ][1] = x2[i] + (h2 / 2);
    //         cout << xi[i * n + j][0]<<"   "<<xi[i * n + j][1] << endl;
    //    }
    //   
    //    
    //}
  /*  for (int i1 = 0; i1 < n; i1++) {
        for (int j1 = 0; j1 < n; j1++) {
            i = i1 * n + j1;
            for (int i2=0; i2 < n; i2++) {
                for (int j2=0; j2 < n; j2++) {
                    j = j2 + i2 * n;
                    A[i][j] = del(i, j) - lymda * middlepryam2(x1[i1], x1[j],x2[i2],x2[j2], xi[i][0], xi[j][1]);
                }

            }
            

        
        }*/

        for (i = 0; i < n * n; i++){ 
            int i1 = i / n;
            int j1 = i % n;

            for ( j = 0; j < n*n; j++){ 
          
                int i2 = j / n;
                int j2 = j % n;

                A[i][j] = del(i,j) - lymda * middlepryam2(x1[i2], x1[i2+1], x2[j2], x2[j2+1], x1[i1]+h1/2, x2[j1] + h2 / 2);

            }
            A[i][n * n] = 1;
            //1 exp(icomp* k0*(x1[i1] + h1 / 2))
        
        }



    //}
    
        /*A[i][n*n] exp(icomp*1000.0*xi[i][0]*xi[i][1]);1;=*/
        
   
    //
    //for (i = 0; i < n*n; i++) {
    //    for (j = 0; j < n*n; j++) {
    //        A[i][j] = del(i, j) - lymda * middlepryam2(x1[i], x1[i+1],x2[j],x2[j+1], xi[i * n + j][0], xi[i * n + j][1]);

    //    }
    //}



    //for (i = 0; i < n2; i++) {
    //    for (j = 0; j < n2 + 1; j++) {

    //        cout << A[i][j] << " ";


    //    }
    //    cout << endl;


    //}

    Gauss(0, A);


    for (i = 0; i < n * n; i++) {
        //for (j = 0; j < n * n + 1 ; j++) {

        //    cout << A[i][j] << " ";

        //}
       cc[i] = A[i][n * n];     
    /*   cout << cc[i] << endl;*/
       /* cout << endl;*/
    }

    /*
    cout << 'ded' << endl;*/

  for (i = 0; i < n; i++) {
      int i1 = j / n;
      int j1 = j % n;
        for (j = 0; j < n; j++) {
       
        out1 << x1[i] << " " << x2[j] << " " << abs(u(x1[i]+h1/2, x2[j]+h2/2, cc, x1, x2)) << endl;
        }
    }

   /*   cout << abs(u(t1[i], cc, t1, t2)) << " " << endl;
    res1[i] = u(t1[i], c, t1, t2);

*/

    return 0;
}

