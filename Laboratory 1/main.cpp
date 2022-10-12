#include <iostream>
#include <cmath>
#include "resources/matrix.h"

using namespace std;

matrix function(matrix x) {

    return -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * M_1_PI, 2) + 0.002 * pow(0.1 * x(), 2));
}

int main() {

    matrix m1 = matrix(0.0);
    matrix m2 = matrix(4.0);
    cout << function(m1);

    double x0, xn, alpha;

    cin >> x0 >> xn >> alpha;
    if (function(m1) == function(m2)) {
        cout<<"x_0 = "<<x0<<", x_n = "<<xn;
    }
    else if (function(m1) > function(m2)) {
        d = -d;
        xn = x0 + d;
        if(function(m1) >= function(m2)){
            cout<<"x_0 = "<<x0<<", x_n - d = "<<xn-d;
        }
    }
    while (function(m1) < function(m2)){

    }

        return 0;
}