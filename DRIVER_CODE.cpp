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

#include "libs/matrix.h"
#include"libs/opt_alg.h"
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
    double* interval=new double[2];
    solution fibSol,lagSol;
    //ekspansja- dobre wyniki inne przedzia³y
//    interval= expansion(func_lab_1,-100.0,1.0,2,1000);
//    printf("[%f,%f]",interval[0],interval[1]);
    //fibonacci- ok
//    fibSol= fib(func_lab_1,10,100,0.00001);
//    cout<<"("<<fibSol.x<<","<<fibSol.y<<")"<<endl;
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