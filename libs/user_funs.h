#pragma once

#include"ode_solver.h"


matrix problem(double t, matrix Y, matrix ud1, matrix ud2);

matrix Fr(matrix X, matrix ud1, matrix ud2);

matrix func_lab_2(matrix x, matrix ud1, matrix ud2);

matrix func_lab_3(matrix x, matrix ud1, matrix ud2);

matrix ff2R(matrix x, matrix ud1, matrix ud2);

matrix df2(double t, matrix Y, matrix ud1, matrix ud2);

matrix df3(double t, matrix Y, matrix ud1, matrix ud2);

matrix Fr3(matrix X, matrix ud1, matrix ud2);

matrix fT3(matrix x, matrix ud1, matrix ud2);

