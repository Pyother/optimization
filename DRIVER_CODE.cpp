/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
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

    //ekspansja- dobre wyniki inne przedzia�y
    interval= expansion(func_lab_1,10,1.0,1.5,1000);
    printf("[%f,%f]",interval[0],interval[1]);

    //fibonacci- ok
//    fibSol= fib(func_lab_1,10,100,0.00001);
//    cout<<"("<<fibSol.x<<","<<fibSol.y<<")"<<endl;

    //Lagrange- ok ale wi�cej wywo�a� funkcji (niepotrzebne wywo�ania?)
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