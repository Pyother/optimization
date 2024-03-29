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
    matrix *Y = solve_ode(problem, 0, 1, 1000, Y0, ud1, X);

    int len = get_len(Y[0]);
    double max = Y[1](0, 2);
    for (int i = 0; i < len; i++)
        if (max < Y[1](i, 2))
            max = Y[1](i, 2);
    double y = abs(max - 50);
    return y;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0(2, 1), Y_ref(2, new double[2]{3.14, 0});
    matrix *Y = solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, x);
    int n = get_len(Y[0]);
    y = 0;
    for (int i = 0; i < n; i++)
        y = y + 10 * pow(Y_ref(0) - Y[1](i, 0), 2) +
            pow(Y_ref(1) - Y[1](i, 1), 2) +
            pow(x(0) * (Y_ref(0) - Y[1](i, 0)) + x(1) * (Y_ref(1) - Y[1](i, 1)), 2);
    y = y * 0.1;
    return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {

    double mr = 1, mc = 9, l = 0.5, b = 0.5;
    double I = mr * l * l / 3 + mc * l * l;
    matrix dY(2, 1);
    dY(0) = Y(1);
    dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;
    return dY;
}


matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {
    double c = 0.47, r = 0.12, m = 0.6, r0 = 1.2, g = 9.81;
    double S = 3.14 * r * r, w = 24.9;
    double Dx = 0.5 * c * r0 * S * Y(1) * abs(Y(1));
    double Dy = 0.5 * c * r0 * S * Y(3) * abs(Y(3));
    double FMx = r0 * Y(3) * w * 3.14 * pow(r, 3);
    double FMy = r0 * Y(1) * w * 3.14 * pow(r, 3);


    matrix dY(4, 1);

    dY(0) = Y(1);
    dY(1) = (-Dx - FMx) / m;
    dY(2) = Y(3);
    dY(3) = (-m * g - Dy - FMy) / m;
    return dY;
}

matrix fT3(matrix x, matrix ud1, matrix ud2) {
    double c = 0, d = 0;
    matrix y = func_lab_3(x, ud1, ud2);
    if (-x(0) + 1 > 0)
        y = 1e10;
    else
        y = y + ud2 * pow(-x(0) + 1, 2);
    if (-x(1) + 1 > 0)
        y = 1e10;
    else
        y = y + ud2 * pow(-x(1) + 1, 2);

    if (norm(x) - a > 0) {
        y = y + ud2 * pow(norm(x) - ud1, 2);
    }
    return y;
}

matrix Fr3(matrix X, matrix ud1, matrix ud2) {

    matrix y;
    matrix Y0 = matrix(4, new double[4]{0, X(0), 100, 0});
    matrix *Y = solve_ode(df3, 0, 0.01, 7, Y0, ud1, X(1));

    int n = get_len(Y[0]);
    int i50 = 0, i0 = 0;
    for (int i = 0; i < n; i++) {
        if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
            i50 = i;
        if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
            i0 = i;
        y = -Y[1](i0, 0);
        if (abs(X(0)) - 10 > 0)
            y = y + ud2 * pow(abs(X(0)) - 10, 2);
        if (abs(X(1)) - 25 > 0)
            y = y + ud2 * pow(abs(X(1)) - 25, 2);
        if (abs(Y[1](i50, 0) - 5) - 1 > 0)
            y = y + pow(abs(Y[1](i50, 0) - 5) - 1, 2);
    }


    return y;
}


//matrix F4T(matrix X, matrix ud1, matrix ud2) {
//    matrix y;
//    if (norm(ud2(0, 0)))
//        y = pow(X(0) + 2 * X * (1));
//    else
//        y = F4T(ud2[0] + X * ud2[1], ud1);
//    return y;
//
//}

matrix gf(matrix X,matrix ud1,matrix ud2){
    matrix g(2,1);
    g(0)=10*X(0)+8*X(1)-34;
    g(1)=8*X(0)+10*X(1)-38;
    return g;
}
matrix Hf(matrix X,matrix ud1,matrix ud2){
    matrix h(2,2);
    h(0,0)=10;
    h(1,0)=8;
    h(0,1)=8;
    h(1,1)=10;
    return h;
}

//matrix FR4(matrix teta,matrix ud1, matrix)


//matrix F4T(matrix X, matrix ud1, matrix ud2) {
//
//    if(norm(ud2(0,0)))
//        y=pow()
//}
matrix fR4(matrix teta, matrix ud1, matrix ud2){
    matrix y;
    int m = 100;
    int n = 3;//get_len(x);
    static matrix X(n, m);
    static matrix Y(1, m);

    //czytamy pliki XData.txt>>X

    double h;
    y = 0;
    for (int i = 0; i < m; i++){
        h = m2d(trans(teta) * X(i));
        h = 1.0 / (1/0 + exp(-h));
        y = y - Y(0, i) * log(h) - (1 - Y(0, i) * log(1-h));
    }
    y = y / m;
    return y;
}

matrix gfR4(matrix teta, matrix ud1, matrix ud2){
    int m = 100;
    int n = get_len(teta);
    static matrix g(n, 1);

    matrix X(100),Y;



    double h;
    for (int j = 0; j < n; j++){
        for (int i = 0; i < m; i++){
            h  = m2d(trans(teta) * X[i]);
            h = 1 / (1 +exp(-h));
            g(j) = g(j) + X(j, i) * (h - Y(0, i));
        }
        g(j) = g(j) / m;
    }
    return g;
}

matrix df6(double t, matrix Y, matrix ud1, matrix ud2) {
    double m1 = 5, m2 = 5, k1 = 1, k2 = 1, F = 1;
    double b1 = ud2(0), b2 = ud2(1);
    matrix dY(4, 1);
    dY(0) = Y(1);
    dY(1) = (-b1 * Y(1) - b2 * (Y(1) - Y(3)) - k1 * Y(0) - k2 * (Y(0) - Y(2))) / m1;
    dY(2) = Y(3);
    dY(3) = (F + b2 * (Y(1) - Y(3)) + k2 * (Y(0) - Y(2))) / m2;
}


matrix ff6R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    int N = 1001;
    static matrix X(N, 2);
    static bool read = true;
    if (read) {
        ifstream S("../polozenia.txt"); // Deleted
        S >> X;
        S.close();
        read = false;
    }
    matrix Y0(4, new double[4]{0, 0, 0, 0});
    matrix *Y = solve_ode(df6, 0, 0.1, 100, Y0, ud1, x[0]);
    y = 0;
    for (int i = 0; i < N; ++i) {
        std::cout << i << " " << Y[1](i, 0) << " " << Y[1](i, 2) << std::endl;
        y = y + abs(X(i, 0) - Y[1](i, 0)) + abs(X(i, 1) - Y[1](i, 2));
    }
    y = y / (2.0 * N);
}




// Deleted lab5


matrix ff6T(matrix x, matrix ud1, matrix ud2){
    matrix y;
    y = pow(x(0) + 2*x(1)-7,2) + pow(2 * x(0) + x(1) - 5,2);
    return  y;
}

