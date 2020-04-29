#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double * xt = new double[200];
double * yt = new double[200];

double f(double x, int num) {
    double x1;
    double x2;
    if(num == 0) {
        x1 = x;
        x2 = 10;
    }
    if(num == 1) {
        x1 = 9;
        x2 = x;

    }
    double f = (pow(x1,2)) + 13*x1 +(pow(x2,2)) + 4*x2 + 37;

	return f;
}

double g(double * x) {
double f = (pow(x[0],2)) + 13*x[0] +(pow(x[1],2)) + 4*x[1] + 37;
return f;
}

double paul(double x1, double Dx, int n);

int main() {
    ofstream fout;
    fout.open("pt.dat");
    double * x = new double[2];
    x[0] = 9;
    x[1] = 10;
    fout << x[0] << " " << x[1] << endl;
    cout << 1 << endl;
    for (int i = 0; i < 2; i++) {
    x[i] = paul(x[i], 0.01, i);
    }
    cout << "Method Pauella: x[0] = " << x[0] << " x[1] =  " << x[1] << " func(x1, x2) = " << g(x) << " In iterations " << xt[199] << endl;
    for(int i = 0; i < 16; i++) {
        if(i == 0 && (fabs(xt[i])>0.0001) && (fabs(yt[i])>0.0001)) {
            fout << xt[0] << " " << yt[0] << endl;
        }
        if(i < 16 && i >= 1 && (fabs(xt[i])>0.0001) && (fabs(yt[i])>0.0001)) {
            fout << xt[i] << " " << yt[i]<< endl;
        }
    }
    fout.close();
    system("python yx.py");
return 0;
}

double paul(double x1, double Dx, int n) {
    double * pt = new double[200];
    pt[199] = 0;
    int iter = 0;
    double e1 = 0.1;
    double e2 = 0.1;
    double x, x2, x3, xmin;
    double fx1, fx2, fx3, fmin; 
    double num, denum;
    x2 = x1 + Dx;
    fx1 = f(x1, n);
    fx2 = f(x2, n);
    if (fx1 > fx2) {
        x3 = x1 + 2*Dx;
    } else {
        x3 = x1 - Dx;
    }
    //ШАГ 4
    for ( ;; ) {
        fx3 = f(x3, n);
        if(fx1 < fx2 && fx1 < fx3) {
            fmin = fx1;
            xmin = x1;
        }
        if(fx2 < fx1 && fx2 < fx3) {
            fmin = fx2;
            xmin = x2;
        }
        if(fx3 < fx1 && fx3 < fx2) {
            fmin = fx3;
            xmin = x3;
        }
        //ШАГ 5
        num = (pow(x2, 2)-pow(x3, 2))*fx1+(pow(x3, 2)-pow(x1, 2))*fx2+(pow(x1, 2)-pow(x2, 2))*fx3;
        denum = (x2-x3)*fx1+(x3-x1)*fx2+(x1-x2)*fx3;
        x = 0.5*(num/denum);
        //ШАГ 6
        pt[iter] = x;
        if(n == 0) {
            xt = pt;
        }
        if(n == 1) {
            yt = pt;
        }
        if (pt[199] < iter) {
                pt[199] = iter;
            }
        if(fabs(xmin-x) < e1 && fabs(fmin-f(x, n)) < e2) {
            return x;
        } else {
            iter++;
            if(x < xmin) {
                xmin = x;
                x1 = xmin;
                x2 = x1 + Dx;
                fx1 = f(x1, n);
                fx2 = f(x2, n);
                if (fx1 > fx2) {
                    x3 = x1 + 2*Dx;
                } else {
                    x3 = x1 - Dx;
                }
            } else {
                x1 = xmin;
                x2 = x1 + Dx;
                fx1 = f(x1, n);
                fx2 = f(x2, n);
                if (fx1 > fx2) {
                    x3 = x1 + 2*Dx;
                } else {
                    x3 = x1 - Dx;
                }
            }
        }
    }
    return x;
}
