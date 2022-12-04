#include"../libs/user_funs.h"

const double a = 0.98;
const double b = 0.63;
const double g = 9.81;

const double Pa = 0.75;
const double Pb = 1;

const double Ta = 90;

const double Fbin = 0.01;
const double Tbin = 10;

const double Db = 36.5665 / 10000;

//assuming

matrix problem(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(3, 1);
    double Fa, Fb;
    Fa = a * b * m2d(ud2) * sqrt(2 * g * Y(0) / Pa);
    Fb = a * b * Db * sqrt(2 * g * Y(1) / Pb);

    if (Y(0) < 0) Fa = 0;
    if (Y(1) < 0) Fb = 0;

    dY(0) = -Fa;
    dY(1) = Fa + Fbin - Fb;
    dY(2) = (Fbin / Y(1)) * (Tbin - Y(2)) + (Fa / Y(1)) * (Ta - Y(2));
//    if (ud1(0)&& (int)(t*10)%10==0)
//    {
//        ud1(0)=ud1(0)+dY(0);
//        ud1(1)=ud1(1)+dY(1);
//        ud1(2)=ud1(2)+dY(2);
//       cout <<"t"<<t<< " A:" << Y(0) + dY(0) << ",B: " << Y(0) + dY(1) << ", temp:" << Y(0) + dY(2) << endl;
//    }
    return dY;
}


matrix Fr(matrix X, matrix ud1, matrix ud2) {
//cout<<ud2(0);
    matrix Y0 = matrix(3, new double[3]{5, 1, 10});
    matrix *Y = solve_ode(problem, 0, 1,1000, Y0, ud1, X);

    int len = get_len(Y[0]);
    double max = Y[1](0, 2);
    for (int i = 0; i < len; i++)
        if (max < Y[1](i, 2))
            max = Y[1](i, 2);
    double y = abs(max - 50);
    return y;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
    matrix y;
    matrix Y0(2, 1), Y_ref(2, new double[2]{ 3.14,0 });
    matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, x);
    int n = get_len(Y[0]);
    y = 0;
    for (int i = 0; i < n; ++i)
        y = y + 10 * pow(Y_ref(0) - Y[1](i, 0), 2) +
            pow(Y_ref(1) - Y[1](i, 1), 2) +
            pow(x(0) * (Y_ref(0) - Y[1](i, 0)) + x(1) * (Y_ref(1) - Y[1](i, 1)), 2);
    y = y * 0.1;
    return y;
}


matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
    double mr = 1, mc = 5, l = 0.5, b = 0.5;
    double I = mr * l * l / 3 + mc * l * l;
    matrix dY(2, 1);
    dY(0) = Y(1);
    dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;
    return dY;
}