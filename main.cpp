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
#include "matrix.h"
#include"opt_alg.h"

void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
    cout << "dzialam1\n";
	try
	{

	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
    cout << "dzialam2\n";
	system("pause");
	return 0;
}

void lab1()
{

}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}

// ##########################################
// LAB1 - additional functions
// ##########################################
matrix function(matrix x) {
    return -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * M_1_PI, 2) + 0.002 * pow(0.1 * x(), 2));
}

void XD(){
//    matrix m1 = matrix(0.0);
//    matrix m2 = matrix(4.0);
//    cout << function(m1);
//
//    double x0, xn, alpha;
//
//    cin >> x0 >> xn >> alpha;
//    if (function(m1) == function(m2)) {
//        cout<<"x_0 = "<<x0<<", x_n = "<<xn;
//    }
//    else if (function(m1) > function(m2)) {
//        d = -d;
//        xn = x0 + d;
//        if(function(m1) >= function(m2)){
//            cout<<"x_0 = "<<x0<<", x_n - d = "<<xn-d;
//        }
//    }
//    while (function(m1) < function(m2)){
//
//    }

}
// ##########################################