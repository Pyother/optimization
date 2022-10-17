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
    cout << "dzialam1\n";
    try {
        lab1();
    }
    catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    cout << "dzialam2\n";
    system("pause");
    return 0;
}

void lab1() {
   cout<< *expansion(func_lab_1,0.0,1.0,2.0,50)<<std::endl;

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
    return -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * M_1_PI, 2) + 0.002 * pow(0.1 * x(), 2));
}
// ##########################################