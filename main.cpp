#define _USE_MATH_DEFINES

#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>
using namespace std;


double Lagrange(valarray<double>* data, double param) {
	double tmp, LagPol = 0;
	unsigned n = data[0].size();
	for (unsigned i = 0; i < n; ++i) {
		tmp = data[1][i];
		for (unsigned j = 0; j < n; ++j) {
			if (i == j) continue;
			tmp *= (param - data[0][j]) / (data[0][i] - data[0][j]);
		}
		LagPol += tmp;
	}
	return LagPol;
}

double Newton(valarray<double>* data, double param) {
	unsigned size = data[0].size(), i;
	vector<vector<double>> f(size - 1);
	double ans, tmp;
	for (i = 0; i < size - 1; ++i) {
		f[0].push_back((data[1][i + 1] - data[1][i]) / (data[0][i + 1] - data[0][i]));
	}
	for (unsigned i = 1; i < size - 1; ++i) {
		for (unsigned j = 0; j < size - i - 1; ++j) {
			f[i].push_back((f[i - 1][j + 1] - f[i - 1][j]) / (data[0][j + i + 1] - data[0][j]));
		}
	}
	ans = data[1][0];
	for (unsigned i = 0; i < size - 1; ++i) {
		tmp = f[i].front();
		for (unsigned j = 0; j <= i; ++j) {
			tmp *= (param - data[0][j]);
		}
		ans += tmp;
	}
	return ans;
}

void NewtonPolinom(valarray<double>* data) {
	unsigned size = data[0].size();
	vector<vector<double>> f(size - 1);
	for (unsigned i = 0; i < size - 1; ++i) {
		f[0].push_back((data[1][i + 1] - data[1][i]) / (data[0][i + 1] - data[0][i]));
	}
	for (unsigned i = 1; i < size - 1; ++i) {
		for (unsigned j = 0; j < size - i - 1; ++j) {
			f[i].push_back((f[i - 1][j + 1] - f[i - 1][j]) / (data[0][j + i + 1] - data[0][j]));
		}
	}
	cout << "N(x) = " << data[1][0];
	for (unsigned i = 0; i < size - 1; ++i) {
		cout << " + (" << f[i][0]<<")";
		for (unsigned j = 0; j <= i; ++j) {
			cout << "*(x - " << data[0][j] << ")";
		}
	}
}

void LagrangePolinom(valarray<double>* data) {
	double tmp, LagPol = 0;
	unsigned n = data[0].size();
    cout << "L(x) = ";
	for (unsigned i = 0; i < n; ++i) {
        if(i>0) cout <<" + ";
        cout << data[1][i];
		for (unsigned j = 0; j < n; ++j) {
			if (i == j) continue;
            cout << "*(x - "<<data[0][j]<<"/("<<data[0][i]<<" - "<<data[0][j]<<")";
		}
	}
}

void Gaus(vector<valarray<double>> &a)
{
    size_t n = a.size();
    double k1, k2;
    //обход алгоритма по всей главной диагонали
    for (size_t i = 0; i < n; ++i)
    {
        //зануляем элементы ниже главной диагонали
        for (size_t j = i + 1; j < n; ++j)
        {
            k1 = a[j][i] / a[i][i];
            a[j] -= a[i] * k1;
        }
        //зануляем элементы выше главной диагонали
        for (size_t j = 0; j < i; ++j)
        {
            k2 = a[j][i] / a[i][i];
            a[j] -= a[i] * k2;
        }
    }
    //сокращаем
    for (size_t i = 0; i < n; ++i)
    {
        a[i][n] /= a[i][i];
        a[i][i] = 1;
    }
}

int main(){
    auto f = [](double x){return tan(x) + x;};
    valarray<double> data[2];
    data[0] = {0, M_PI/8, M_PI/4, 3*M_PI/8};
    data[1] = data[0];
    for(auto& el : data[1]) el = f(el);
    NewtonPolinom(data);
    cout<<endl;
    LagrangePolinom(data);
    size_t n = data[0].size();
    vector<valarray<double>> matrix(n);
    for(size_t i = 0; i < n; ++i){
        matrix[i].resize(n + 1);
        for(size_t j = 0; j < n; ++j){
            matrix[i][j] = pow(data[0][i], j);
        }
        matrix[i][n] = data[1][i];
    }
    Gaus(matrix);
    cout<<"K(x) = ";
    for( size_t i = 0; i < n; ++i){
        if(i>0) cout <<" + ";
        cout<<matrix[i][n]<<"*(x^"<<i<<")";
    }
    cout<<endl;
    return 0;
}