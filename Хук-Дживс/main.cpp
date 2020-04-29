#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <fstream>
#include <cstdlib>

using namespace std;
double PI = 3.1415;

double norm(double *x, int n = 2) {
	double r = 0;
	for (int i = 0; i < n; i++)
		r += x[i] * x[i];
	return sqrt(r);
}

double dot(double *a, double *b, int n = 2) {
	double r = 0;
	for (int i = 0; i < n; i++)
		r += a[i] * b[i];
	return r;
}

double f(double *x) {
	return (x[0]*x[0]+13*x[0]+x[1]*x[1]+4*x[1]+37);
}

double dfdx(double *x) {
	return 2 * x[0] + 13;
}

double dfdy(double *x) {
	return 2 * x[1] + 4;
}

double f_s(double *x, double a, double *p) {
	double x_new[] = { x[0] + a*p[0], x[1] + a*p[1] };
	return f(x_new);
}

double svn(double *x, double *p) { // Алгоритм Свенна
	double a = 0;
	double h = 0.05;
	while (1) {
		if (f_s(x, a - h, p) >= f_s(x, a, p) && f_s(x, a, p) <= f_s(x, a + h, p)) {
			break;
		}
		else if (f_s(x, a - h, p) >= f_s(x, a, p) && f_s(x, a, p) >= f_s(x, a + h, p)) {
			a += h / 2.;
		}
		else if (f_s(x, a - h, p) < f_s(x, a, p) && f_s(x, a, p) < f_s(x, a + h, p)) {
			a -= h / 2.;
		}
		else {
			a += h;
		}
	}
	return a;
}

double mki(double *x, double *p) { // Метод квадратичной интерполяции
	double h = 0.025;
	double a = svn(x, p);
	double a1 = a, a2, a3, a_min, a_star;
	do {
		a2 = a1 + h;
		if (f_s(x, a1, p) > f_s(x, a2, p))
			a3 = a1 + 2 * h;
		else
			a3 = a1 - h;

		a_star = 0.5*((a2*a2 - a3*a3)*f_s(x, a1, p) + (a3*a3 - a1*a1)*f_s(x, a2, p) + (a1*a1 - a2*a2)*f_s(x, a3, p));
		double temp = (a2 - a3)*f_s(x, a1, p) + (a3 - a1)*f_s(x, a2, p) + (a1 - a2)*f_s(x, a3, p);

		if (f_s(x, a1, p) <= f_s(x, a2, p) && f_s(x, a2, p) <= f_s(x, a3, p))
			a_min = a1;
		else if (f_s(x, a2, p) <= f_s(x, a1, p) && f_s(x, a1, p) <= f_s(x, a3, p))
			a_min = a2;
		else
			a_min = a3;

		if (temp == 0) {
			a1 = a_min;
			continue;
		}
		a_star /= temp;

		if (a1 <= a_star && a_star <= a3) {
			if (f_s(x, a_min, p) < f_s(x, a_star, p))
				a1 = a_min;
			else
				a1 = a_star;
		}
		else
			a1 = a_star;
	} while (abs(a_min - a_star) >= 1e-4);
	return a_star;
}


double f_h(double *x, double h, int p) {
	double x_new[] = { x[0], x[1] };
	if (p == 0)
		x_new[0] += h;
	else if (p == 1)
		x_new[1] += h;
	else
		x_new[2] += h;
	return f(x_new);
}

void hooke_jeeves(double *x0, double h = 1, double l = 0.9, double eps = 1e-5) {
	ofstream traekt("hj.txt", ios_base::out);
	double x[] = { x0[0], x0[1] };
	double xn[] = { x[0], x[1] };
	int i = 0;
	while (h > eps) {
		xn[0] = x[0]; xn[1] = x[1];
		if (f_h(xn, h, 0) < f(xn))
			xn[0] += h;
		else if (f_h(xn, -h, 0) < f(xn))
			xn[0] -= h;

		if (f_h(xn, h, 1) < f(xn))
			xn[1] += h;
		else if (f_h(xn, -h, 1) < f(xn))
			xn[1] -= h;

		if (xn[0] - x[0] == 0 && xn[1] - x[1] == 0 ) {
			h *= l;
			continue;
		}
		traekt << x[0] << " " << x[1] << endl;
		i++;
		double p[] = { xn[0] - x[0], xn[1] - x[1], 0 };
		double t = norm(p);
		for (int i = 0; i < 2; i++) p[i] /= t;
		x[0] += mki(x, p)*p[0];
		x[1] += mki(x, p)*p[1];
	}
	cout << "iter: " << i << endl;
	cout << "xmax: (" << x[0] << " " << x[1] << " "<< ")" << endl;
	cout << "f(" << x[0] << " " << x[1] << " " << ") = " << f(x) << endl;
}




int main()
{
	double x0[] = { 9., 10. };
	cout << "Hooke-Jeeves method:\n";
	hooke_jeeves(x0);
	system("pause");
	system("python yx.py");
	return 0;
}
