#pragma once

#include"ode_solver.h"


matrix fun(double t, matrix Y, matrix ud1, matrix ud2);

matrix Fr(matrix X, matrix ud1, matrix ud2);