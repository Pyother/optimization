#include"../libs/user_funs.h"


matrix fun(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY = matrix(3, 1);

    dY[0] = a * b * m2d(ud1) * sqrt(2 * g * Y[0] / Pa);
    dY[1] = a * b * m2d(ud1) * sqrt(2 * g * Y[1] / Pb)

    dY[2] =( Fbin/dY[0]) *(Tbin-Y[2])+(Fa / Y[1]) *(Tbin-Y[2])
}


matrix Fr(matrix X, matrix ud1, matrix ud2) {

    matrix Y0 = matrix(3, new double[3]{5, 1, 10});
    matrix *Y = solve_ode(NULL, 0, 1, 1000, ud1, ud2);


}
