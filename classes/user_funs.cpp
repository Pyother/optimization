#include"../libs/user_funs.h"

const double a = 0.98;
const double b = 0.63;
const double g = 9.81;

const double Pa = 0.75;
const double Pb = 1;

const double Ta = 90;

const double Fbin = 0.01;
const double Tbin = 10;

const double Db= 36.5665;

//assuming

matrix fun(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY = matrix(3, 1);

// assuming that Y=volume
    auto Fa = a * b * m2d(ud1) * sqrt(2 * g * m2d(Y[0]) / Pa);
// assuming Db=ud2
    auto Fb = a * b * Db * sqrt(2 * g * m2d(Y[1]) / Pb);

    dY[0] = -Fa;//=Va
    dY[1] = Fa + Fbin - Fb;//=Vb
// there was dY[0] in first parenthesis before
    dY[2] = (Fbin / Y[1]) * (Tbin - Y[2]) + (Fa / Y[1]) * (Ta - Y[2]); //=Tb
//fun(t+1,dY,ud1,ud2);
    return dY;
}


matrix Fr(matrix X, matrix ud1, matrix ud2) {

    matrix Y0 = matrix(3, new double[3]{5, 1, 10});
    matrix *Y = solve_ode(NULL, 0, 1, 1000, ud1, X);

    //double y= abs()

    //return y;
}
