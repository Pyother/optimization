/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
*********************************************/
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <cstdlib>
#include <fstream>
#include "libs/matrix.h"
#include "libs/opt_alg.h"
#include "libs/solution.h"



using namespace std;

matrix func_lab_1(matrix x, matrix ud1, matrix ud2);

void lab1();

void lab2();

void lab3();

void lab4();

void lab5();

void lab6();

int main() {
    try {
        lab1();
    }
    catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    return 0;
}

void lab1() {

    #define tab1 1
    #define tab2 1
    double epsilon = 1e-5;
    double gamma = 1e-200;

    double* interval;
    solution fibSol, lagSol;
    double* Y = new double[3];
    Y[0] = 1;
    Y[1] = 1;
    Y[2] = 1;
    matrix ud2;
    //ekspansja- dobre wyniki inne przedzia³y
    int x0, d = 1, alpha = 2, Nmax = 1000;


    #if tab1

        ofstream expansionFile;
        expansionFile.open("out.csv", ofstream::out);

        for (int i = 0; i < 100; i++) {
            x0 = rand() % 200 + 1;
            interval = expansion(func_lab_1, x0, d, alpha, Nmax, matrix(3, Y), ud2);
            expansionFile << x0 << ", " << interval[0] << ", " << interval[1] << ", " << solution::f_calls << endl;
            solution::clear_calls();
            fibSol = fib(func_lab_1,interval[0],interval[1],0.00001);
            solution::clear_calls();
            lagSol = lag(func_lab_1,interval[0],interval[1],0.0001,0.0000001,1000);
            solution::clear_calls();
        }
        expansionFile.close();

    #endif

    #if tab2

        matrix ab_F(1, 1, 200);
        solution opt_f = fib(func_lab_1, -100, 100, epsilon);
        cout << endl << endl;
        cout << "Fibonacci:" << endl;
        cout << opt_f ;
        cout << "ab_F = " << endl << ab_F << endl;
        cout << endl << endl;

        matrix ab_L(1, 1, 200);
        solution opt_l = lag(func_lab_1,-100,100, epsilon, gamma,1000);
        cout << "Lagrange:" << endl;
        cout << opt_l;
        cout << "ab_F = " << endl << ab_L << endl;

    #endif


    //printf("[%f,%f]",interval[0],interval[1]);

    //fibonacci- ok
//    fibSol= fib(func_lab_1,10,100,0.00001);
//    cout << fibSol.x << endl << fibSol.y << endl;


    //Lagrange- ok ale wiêcej wywo³añ funkcji (niepotrzebne wywo³ania?)
//    lagSol=lag(func_lab_1,-10,1,0.0001,0.0000001,1000);
//    cout<<"("<<lagSol.x<<","<<lagSol.y<<") "<<"calls: "<<lagSol.f_calls<<" flag: "<<lagSol.flag<<endl;
}

void lab2() {

}

void lab3() {

}

void lab4() {

}

void lab5() {

}

void lab6() {

}

// ##########################################
// LAB1 - mathematical functions
// ##########################################
//
//
// PLIK opt_algo.cpp
//
//
// ##########################################
matrix func_lab_1(matrix x, matrix ud1, matrix ud2) {
    return -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
}
// ##########################################