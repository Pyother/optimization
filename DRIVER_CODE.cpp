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
#include "libs/user_funs.h"


using namespace std;

matrix func_lab_1(matrix x, matrix ud1, matrix ud2);

matrix func_lab_2(matrix x, matrix ud1, matrix ud2);

matrix func_lab_3(matrix x, matrix ud1, matrix ud2);
matrix func_lab_3_test(matrix x, matrix ud1, matrix ud2);


void simulation(matrix Da, matrix ud1, matrix ud2);

void lab1();

void lab2();

void lab3();

void lab4();

void lab5();

void lab6();


int main() {
    try {
//        lab1();
//        lab2();
        lab3();
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
    matrix ud1, ud2;

    double Da0 = 1.0e-4;
    double d = 0.01;
    double alpha = 2;
    int Nmax = 1000;
    solution fibSol, lagSol;
//    //ud1 = matrix(3, new double[3]{5, 1, 10});
//    double *Range = expansion(Fr, Da0, d, alpha, Nmax, ud1);
//    cout << "test range: (" << Range[0] << "," << Range[1] << ")" << endl;


//    solution::clear_calls();
//    fibSol = fib(Fr, Range[0], Range[1], epsilon);
//    cout << fibSol.x << " " << fibSol.y << " " << solution::f_calls << endl;
//    fibSol.clear_calls();
//    lagSol = lag(Fr, Range[0], Range[1], epsilon, gamma, 1000);
//    cout << lagSol.x << " " << lagSol.y << " " << solution::f_calls << endl;
//
//     ofstream tabela3File;
//    tabela3File.open("tabela3File.csv", ofstream::out);
//    solution::clear_calls();
//    matrix Da = matrix(1, new double[1]{0.00208884});
//    matrix Y0 = matrix(3, new double[3]{5, 1, 10});
//    matrix *Y = solve_ode(problem, 0, 1, 300, Y0, 0, lagSol.x);
//    tabela3File << Y[0] << " " <<Y[1][2] << " " << endl;
//
//
//    tabela3File.close();
//    //simulation(Da, ud1, ud1);

////    solution fibSol,lagSol;
////    //ekspansja- dobre wyniki inne przedzia?y
//    interval= expansion(func_lab_1,-100.0,1.0,2,1000);
//    printf("[%f,%f]",interval[0],interval[1]);
//    //double* interval;
// solution fibSol, lagSol;
    double *Y = new double[3];
    Y[0] = 1;
    Y[1] = 1;
    Y[2] = 1;
//    matrix ud2;
    //ekspansja- dobre wyniki inne przedzia³y
    int x0;
    d = 1, alpha = 3.342, Nmax = 1000;

    double *interval = new double[2];
#if tab1

    ofstream expansionFile;
    ofstream fibbonacciFile, lagrangeFile;
    expansionFile.open("expansion1.csv", ofstream::out);
    fibbonacciFile.open("fibbonacci1.csv", ofstream::out);
    lagrangeFile.open("lagrange1.csv", ofstream::out);

    for (int i = 0; i < 100; i++) {
        x0 = rand() % 200 + 1;
        interval = expansion(func_lab_1, x0, d, alpha, Nmax, matrix(3, Y), ud2);
        expansionFile << x0 << ", " << interval[0] << ", " << interval[1] << ", " << solution::f_calls << endl;
        solution::clear_calls();
        fibSol = fib(func_lab_1, interval[0], interval[1], 0.00001);
        fibbonacciFile << fibSol.x << " " << fibSol.y << " " << solution::f_calls << endl;
        solution::clear_calls();
        lagSol = lag(func_lab_1, interval[0], interval[1], 0.0001, 0.0000001, 1000);
        lagrangeFile << lagSol.x << " " << lagSol.y << " " << solution::f_calls << endl;
        solution::clear_calls();
    }
    expansionFile.close();
    fibbonacciFile.close();
    lagrangeFile.close();

#endif

#if tab2

    matrix ab_F(1, 1, 200);
    solution opt_f = fib(func_lab_1, interval[0], interval[1], epsilon);
    cout << endl << endl;
    cout << "Fibonacci:" << endl;
    cout << opt_f;
    cout << "ab_F = " << endl << ab_F << endl;
    cout << endl << endl;

    matrix ab_L(1, 1, 200);
    solution opt_l = lag(func_lab_1, interval[0], interval[1], epsilon, gamma, 1000);
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
    matrix ud1 = matrix(2,new double [2]{3.14,0}), ud2 = matrix(2,new double [2]{1,1});

    double s = 0.2, epsilon = 1e-5, alphaHJ = 0.5, alphaRosen = 2, beta = 0.5;
    int Nmax = 1000;
    matrix x0(2, new double[2]{0, 0});
    matrix s_rosen(2, 1, s);
    cout << HJ(func_lab_2, x0, s, alphaHJ, epsilon, Nmax, ud1, ud2) << endl;
    solution::clear_calls();
//    cout << RosenFileTab1(func_lab_2, x0, s_rosen, alphaRosen, beta, epsilon, Nmax, ud1, ud2) << endl;

//    printf("TABELA 1");
//    solution intervalHJ, intervalRosen;
//    ofstream HooksJeevesFileTab1;
//    ofstream RosenFileTab1;
//    HooksJeevesFileTab1.open("HooksJeevesFileTab1.csv", ofstream::out);
//    RosenFileTab1.open("RosenFileTab1.csv", ofstream::out);
//    solution HJSol, RosSol;
//    for (int i = 0; i < 100; i++) {
//        std::default_random_engine rnd1{std::random_device{}()};
//        std::uniform_real_distribution<double> dist(-1, 1);
//        double r1 = dist(rnd1);
//        std::default_random_engine rnd2{std::random_device{}()};
//        std::uniform_real_distribution<double> dist2(-1, 1);
//        double r2 = dist(rnd2);
//        x0 = matrix(2, new double [2] {r1, r2});
//
//        intervalHJ = HJ(func_lab_2, x0, s, alphaHJ, epsilon, Nmax, ud1, ud2);
//        HooksJeevesFileTab1 << x0(0) << ", "  << x0(1) << ", " << intervalHJ.x(0) << ", " << intervalHJ.x(1) << ", "
//                << intervalHJ.y(0) << ", " << solution::f_calls << endl;
//        solution::clear_calls();
//
//        intervalRosen = Rosen(func_lab_2, x0, s_rosen, alphaRosen, beta, epsilon, Nmax, ud1, ud2);
//        RosenFileTab1 << intervalRosen.x(0) << ", " << intervalRosen.x(1) << ", "
//                            << intervalRosen.y(0) << ", " << solution::f_calls << endl;
//        solution::clear_calls();
//    }
//    HooksJeevesFileTab1.close();
//    RosenFileTab1.close();


    printf("\n\n\n");
    solution rosenSol=Rosen(func_lab_2, x0, s_rosen, alphaRosen, beta, epsilon, Nmax, ud1, ud2);
    cout << rosenSol;
    solution::clear_calls();
    solution hjSol=HJ(func_lab_2, x0, s, alphaHJ, epsilon, Nmax, ud1, ud2);
    cout << hjSol;


    matrix Y0(2, 1), Y_ref(2, new double[2]{ 3.14,0 });
    matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, rosenSol.x);
    ofstream HJFileTabSim;
    HJFileTabSim.open("HJFileTabSim.csv", ofstream::out);
    HJFileTabSim<<Y[1]<<endl;
    HJFileTabSim.close();
}

void lab3() {
    solution simp;
    matrix x=matrix(2,new double[2]{-3,2}),ud1,ud2;
    double s=3,aplha=1.0,beta=0.5,gamma=0.5,delta=0.5;

    
   // matrix ud1 = matrix(2,new double [2]{3.14,0}), ud2 = matrix(2,new double [2]{1,1});
    matrix x0zew(2, new double[2]{1, 2});
    ofstream S("lab_3.csv");
    double c = 1, dc = 0.2, epsilon = 1e-3;
    int Nmax = 10000;
    matrix a = 5;
    S << "X0" << ";" << " " << "X1" << std::endl << std::endl;

    S << "Sym_MN dc 0.2 " << endl << endl;
    S << "x1*" << ";" << " " << "x2*" << ";" << ";" << " " << "y*" << ";" <<
      " " << "f_calls*" << ";" << " " << "r" << std::endl;

    // zewnetrzena dla dc > 1
    // wewnetrzna dla dc < 1
    for (int i = 0; i < 100; i++)
    {
        solution test = pen(fT3, x0zew, c, dc, epsilon, Nmax, ud1, ud2);
        S << x0zew(0) << ", " << x0zew(1) << ", " << test.x(0) << ", " << test.x(1) << ", "
                << test.y(0) << ", " << solution::f_calls << endl;
        solution::clear_calls();
    }

    dc = 2;
    S << "Sym_MN dc 2" << endl << endl;
    S << "x1*" << ";" << " " << "x2*" << ";" << ";" << " " << "y*" << ";" <<
      " " << "f_calls*" << ";" << " " << "r" << std::endl;
    for (int i = 0; i < 100; i++)
    {
        solution test = pen(fT3, x0zew, c, dc, epsilon, Nmax, ud1, ud2);
        S << test;
        S << " " << sqrt(pow(test.x(0), 2) + pow(test.x(1), 2)) << endl;
        solution::clear_calls();
    }
    S.close();

    // #################################################################################################################
    // #################################################################################################################
    // #################################################################################################################
    // #################################################################################################################
    // #################################################################################################################
    printf("TABELA 1");
    ofstream zew_penalty, wew_penalty;
    zew_penalty.open("zew_penalty.csv", ofstream::out);
    wew_penalty.open("wew_penalty.csv", ofstream::out);

    ud1 = 4;
    double dcz = 2, dcw = 0.2;

    zew_penalty << "x0zew(0)" << ", " << "x0zew(1)" << ", " << "test.x(0)" << ", " << "test.x(1)" << ", "
                << "test.y(0)" << ", " << "solution::f_calls" << endl;
    wew_penalty << "x0zew(0)" << ", " << "x0zew(1)" << ", " << "test.x(0)" << ", " << "test.x(1)" << ", "
                << "test.y(0)" << ", " << "solution::f_calls" << endl;
    matrix x0wew(2, new double[2]{1, 2});
    for (int i = 0; i < 100; i++){
        std::default_random_engine random1{std::random_device{}()};
        std::uniform_real_distribution<double> v1(-10, 10);
        double val1 = v1(random1);
        std::default_random_engine random2{std::random_device{}()};
        std::uniform_real_distribution<double> v2(-10, 10);
        double val2 = v2(random2);
        x0zew = matrix(2, new double [2] {val1, val2});
        solution pen_func_z = pen(fT3, x0zew, c, dcz, epsilon, Nmax, ud1, ud2);
//        zew_penalty << pen_func_z;
        zew_penalty << "a,  " << x0zew(0) << ", " << x0zew(1) << ", " << pen_func_z.x(0) << ", " << pen_func_z.x(1) << ", "
                    << pen_func_z.y(0) << ", " << solution::f_calls << endl;
        solution::clear_calls();

        x0wew = matrix(2, new double [2] {val1, val2});
        solution pen_func_w = pen(fT3, x0wew, c, dcw, epsilon, Nmax, ud1, ud2);
//        wew_penalty << pen_func_z;
        wew_penalty << "b,  " << x0wew(0) << ", " << x0wew(1) << ", " << pen_func_w.x(0) << ", " << pen_func_w.x(1) << ", "
                    << pen_func_w.y(0) << ", " << solution::f_calls << endl;
        solution::clear_calls();
    }
    zew_penalty.close();
    wew_penalty.close();
}

void lab4() {

}

void lab5() {

}

void lab6() {

}

// ##########################################
// ##########################################
// Mathematical functions for LABs
// ##########################################
// ##########################################
matrix func_lab_1(matrix x, matrix ud1, matrix ud2) {
    return -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
}

matrix func_lab_2(matrix x, matrix ud1, matrix ud2) {
    return (pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2);
}

matrix func_lab_3(matrix x, matrix ud1, matrix ud2) {
    return (sin(3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)))) /
           (3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2)));
}
matrix func_lab_3_test(matrix x, matrix ud1, matrix ud2) {
    return pow(x(0)-3,2)+pow(x(1)-3,2);
}


void simulation(matrix Da, matrix ud1, matrix ud2) {
//cout<<ud2(0);
    solution fibSol, lagSol;

    double epsilon = 1e-5;
    double gamma = 1e-200;
//    matrix Y0 = matrix(3, new double[3]{5, 1, 10});
//    matrix Y = problem(0, Y0, ud1, Da);
//    Y0(0) += Y(0);
//    Y0(1) += Y(1);
//    Y0(2) += Y(2);
//    matrix dY0;
//    for (int t = 0; t < 30; t++) {
//        dY0 = problem(t, Y, ud1, Da);
//        Y0(0) += dY0(0);
//        Y0(1) += dY0(1);
//        Y0(2) += dY0(2);
//        cout << "A:" << Y0(0) << ",B: " << Y0(1) << ", temp:" << Y0(2) << endl;
//    }
    ud1 = matrix(3, new double[3]{5, 1, 10});
    fibSol = fib(Fr, Da(0), Da(0), epsilon);
    cout << fibSol.x << " " << fibSol.y << " " << solution::f_calls << endl;
    solution::clear_calls();
    fibSol = lag(Fr, Da(0), Da(0), epsilon, gamma, 1000, ud1, ud1);
    cout << fibSol.x << " " << fibSol.y << " " << solution::f_calls << endl;




}

// ##########################################
// ##########################################
