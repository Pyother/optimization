#include"../libs/opt_alg.h"
#include"../libs/matrix.h"

double *
expansion(matrix(ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2) {
    auto *p = new double[2]{0, 0};
    solution X0(x0), X1(x0 + d);
    X0.fit_fun(ff, ud1, ud2);
    X1.fit_fun(ff, ud1, ud2);

    if (X0.y == X1.y) {
        p[0] = m2d(X0.x);
        p[1] = m2d(X1.x);
//            cout << "if 1\n";
        cout << " Liczba wywolan  = " << solution::f_calls << endl;
        return p;

    }
    if (X0.y < X1.y) {
        d = -d;
        X1.x = X0.x + d;
        X1.fit_fun(ff, ud1, ud2);
        if (X1.y >= X0.y) {
            p[0] = X1.x();
            p[1] = X0.x() - d;
            cout << " Liczba wywolan  = " << solution::f_calls;
            return p;
        }
    }
    solution X2(X1.x);
    int i = 0;
    while (true) {
        i++;
        X2.x = x0 + pow(alpha, i) * d;
        X2.fit_fun(ff, ud1, ud2);

        if (X2.y >= X1.y || solution::f_calls > Nmax)
            break;
        X0 = X1;
        X1 = X2;

    }
    d > 0 ? p[0] = X0.x(), p[1] = X2.x() : (p[0] = X2.x(), p[1] = X0.x());
    return p;
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        Xopt.ud = b - a;
        int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
        int *F = new int[n]{1, 1};
        for (int i = 2; i < n; ++i)
            F[i] = F[i - 2] + F[i - 1];
        solution A(a), B(b), C, D;
        C.x = B.x - 1.0 * F[n - 2] / F[n - 1] * (B.x - A.x);
        D.x = A.x + B.x - C.x;
        C.fit_fun(ff, ud1, ud2);
        D.fit_fun(ff, ud1, ud2);


        for (int i = 0; i <= n - 3; ++i) {
            if (C.y < D.y) {
                A.x = A.x;
                B.x = D.x;
            } else {
                B.x = B.x;
                A.x = C.x;
            }
            C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
            D.x = A.x + B.x - C.x;

            C.fit_fun(ff, ud1, ud2);
            D.fit_fun(ff, ud1, ud2);
        }
        Xopt = C;
        Xopt.flag = 0;
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }

}

solution
lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;

        Xopt.ud = b - a;
        solution A(a), B(b), C, D, D_old(a);
        C.x = (a + b) / 2;
        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);
        double l, m;
        while (true) {
            l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y * (pow(C.x) - pow(A.x)) + C.y * (pow(A.x) - pow(B.x)));
            m = m2d(A.y * (B.x - C.x) + B.y * (C.x - A.x) + C.y * (A.x - B.x));
            if (m <= 0) {
                Xopt = D_old;
                Xopt.flag = 2;
                return Xopt;
            }
            D.x = 0.5 * l / m;

            D.fit_fun(ff, ud1, ud2);
            if (A.x <= D.x && D.x <= C.x) {
                if (D.y < C.y) {
                    A.x = A.x;
                    B.x = C.x;
                    C.x = D.x;
                } else {
                    A.x = D.x;
                    C.x = C.x;
                    B.x = B.x;
                }

            } else if (C.x <= D.x && D.x <= B.x) {
                if (D.y < C.y) {
                    A.x = C.x;
                    C.x = D.x;
                    B.x = B.x;
                } else {
                    A.x = A.x;
                    C.x = C.x;
                    B.x = D.x;
                }
            } else {
                Xopt = D_old;
                Xopt.flag = 2;
                return Xopt;
            }
            A.fit_fun(ff, ud1, ud2);
            B.fit_fun(ff, ud1, ud2);
            C.fit_fun(ff, ud1, ud2);

            Xopt.ud.add_row((B.x - A.x)());
            if (B.x - A.x < epsilon || abs(D.x() - D_old.x()) < gamma) {
                Xopt = D;
                Xopt.flag = 0;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt = D;
                Xopt.flag = 1;
                break;
            }
            D_old = D;
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution lag(...):\n" + ex_info);
    }
}


solution
HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1,
   matrix ud2) {
    try {
        ofstream HJFileTab3;
        HJFileTab3.open("HJFileTab3.csv", ofstream::out);
        solution XB, XB_old, X;
        XB.x = x0;
        XB.flag = 0;
        XB.fit_fun(ff, ud1, ud2);
        while (true) {

            X = HJ_trial(ff, XB, s, ud1, ud2);
            //cout << X.x(0) << " " << X.x(1) << endl;
            if (X.y < XB.y) {
                while (true) {
                    XB_old = XB;
                    XB = X;
                    HJFileTab3 << X.x(0) << ", " << X.x(1) << endl;
                    X.x = XB.x + XB.x - XB_old.x;

                    X.fit_fun(ff, ud1, ud2);
                    X = HJ_trial(ff, X, s, ud1, ud2);
                    if (X.y >= XB.y)
                        break;
                    if (XB.f_calls > Nmax) {
                        XB.flag = 1;
                        return XB;
                    }
                }
            } else
                s *= alpha;
            if (s < epsilon) {
                XB.flag = 0;
                HJFileTab3.close();
                return XB;
            } else if (XB.f_calls > Nmax) {
                XB.flag = 1;
                HJFileTab3.close();
                return XB;
            }
        }

    }
    catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        int n = get_dim(XB);
        solution X;
        X.flag = 0;
        matrix d = ident_mat(n);
        for (int i = 0; i < n; i++) {
            X.x = XB.x + s * d[i];
            X.fit_fun(ff, ud1, ud2);
            if (X.y < XB.y) {
                XB = X;
            } else {
                X.x = XB.x - s * d[i];
                X.fit_fun(ff, ud1, ud2);
                if (X.y < XB.y) {
                    XB = X;
                }

            }

        }
        return XB;
    }
    catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution
Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax,
      matrix ud1, matrix ud2) {
    try {
        ofstream RosenFileTab3;
        RosenFileTab3.open("RosenFileTab3.csv", ofstream::out);
        solution X(x0), Xt;
        int n = get_dim(X);
        matrix l(n, 1),
                p(n, 1),
                s(s0),
                D = ident_mat(n);
        X.fit_fun(ff, ud1, ud2);
        while (true) {
            for (int i = 0; i < n; ++i) {
                Xt.x = X.x + s(i) * D[i];
                Xt.fit_fun(ff, ud1, ud2);
                RosenFileTab3 << X.x(0) << ", " << X.x(1) << endl;
                if (Xt.y < X.y) {
                    X = Xt;
                    l(i) += s(i);
                    s(i) *= alpha;
                } else {
                    ++p(i);
                    s(i) *= -beta;
                }
            }
            bool change = true;
            for (int i = 0; i < n; ++i)
                if (p(i) == 0 || l(i) == 0) {
                    change = false;
                    break;
                }
            if (change) {
                matrix Q(n, n), v(n, 1);
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j <= i; ++j)
                        Q(i, j) = l(i);
                Q = D * Q;
                v = Q[0] / norm(Q[0]);
                D.set_col(v, 0);
                for (int i = 1; i < n; ++i) {
                    matrix temp(n, 1);
                    for (int j = 0; j < i; ++j)
                        temp = temp + trans(Q[i]) * D[j] * D[j];
                    v = (Q[i] - temp) / norm(Q[i] - temp);
                    D.set_col(v, i);
                }
                s = s0;
                l = matrix(n, 1);
                p = matrix(n, 1);
            }
            double max_s = abs(s(0));
            for (int i = 1; i < n; ++i)
                if (max_s < abs(s(i)))
                    max_s = abs(s(i));
            if (max_s < epsilon) {
                X.flag = 0;
                RosenFileTab3.close();
                return X;
            } else if (solution::f_calls > Nmax) {
                X.flag = 1;
                RosenFileTab3.close();

                return X;
            }
        }
    }
    catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2) {
    double alpha = 4, beta = 0.5, gamma = 2, delta = 0.5, s = 0.5;
    solution X(x0), X1;
    matrix c0(2, new double[2]{c, dc});
    while (true) {
        X1 = sym_NM(ff, X.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c0(0));//ud2=c0?
        if (solution::f_calls > Nmax || norm(X.x - X1.x) < epsilon)
            return X1;
        c0(0) *= dc;
        X = X1;
    }
}

solution
sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta,
       double epsilon, int Nmax, matrix ud1, matrix ud2) {
    int n = get_len(x0);
    matrix D = ident_mat(n);
    int N = n + 1;
    solution *S = new solution[N];
    S[0].x = x0;
    S[0].fit_fun(ff, ud1, ud2);
    for (int i = 1; i < N; ++i) {
        S[i].x = S[0].x + s * D[i - 1];
        S[i].fit_fun(ff, ud1, ud2);
    }
    solution PR, PE, PN;
    matrix pc;
    int i_min, i_max;
    while (true) {
        i_min = i_max = 0;
        for (int i = 1; i < N; ++i) {
            if (S[i_min].y > S[i].y)
                i_min = i;
            if (S[i_max].y < S[i].y)
                i_max = i;
        }
        pc = matrix(n, 1);
        for (int i = 0; i < N; ++i)
            if (i != i_max)
                pc = pc + S[i].x;
        pc = pc / n;
        PR.x = pc + alpha * (pc - S[i_max].x);
        PR.fit_fun(ff, ud1, ud2);
        if (S[i_min].y <= PR.y && PR.y < S[i_max].y)
            S[i_max] = PR;
        else if (PR.y < S[i_min].y) {
            PE.x = pc + gamma * (PR.x - pc);
            PE.fit_fun(ff, ud1, ud2);
            if (PE.y < PR.y)
                S[i_max] = PE;
            else
                S[i_max] = PR;
        } else {
            PN.x = pc + beta * (S[i_max].x - pc);
            PN.fit_fun(ff, ud1, ud2);
            if (PN.y < S[i_max].y)
                S[i_max] = PN;
            else {
                for (int i = 0; i < N; ++i)
                    if (i != i_min) {
                        S[i].x = delta * (S[i].x + S[i_min].x);
                        S[i].fit_fun(ff, ud1, ud2);
                    }
            }
        }
        double max_s = norm(S[i_min].x - S[0].x);
        for (int i = 1; i < N; ++i)
            if (max_s < norm(S[i_min].x - S[i].x))
                max_s = norm(S[i_min].x - S[i].x);
        if (max_s < epsilon || solution::f_calls > Nmax)
            return S[i_min];
    }
}

// #####################################################################################################################
solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1 , matrix ud2 ) {
    try {
        solution Xopt;
        Xopt.ud = trans(x0);
        int n = get_len(x0); //rozmiar problemu
        solution X0, X1;
        X0.x = x0; // pkt startowy
        matrix d(n, 1); // kierunek
        matrix P(n, 2);
        solution h; //dl kroku do wyznaczenia
        double *ab; // zwracane przez funkcje ekspansji
        while (true) {
//        X0.grad(gf, ud1, ud2); // gradient
//        d = -X0.g; //solution.g - gradient, .H - hesjan
            d = -X0.grad(gf, ud1, ud2);
            // czy stalo/zmienno-krokowa
            if (h0 < 0) {
                // met zmiennokrokowa
                P.set_col(X0.x, 0);
                P.set_col(d, 1);
//            P[0] = X0.x;
//            P[1] = d;
//            ab = expansion(ff,0, 1, 1.2, Nmax);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
//            h = golden(ff,ab[0], ab[1], epsilon, Nmax);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
                X1.x = X0.x + h.x * d;
            } else { // stalokrokowa
                X1.x = X0.x + h0 * d;
            }

//        ud->add_row(trans(X1.x));
            Xopt.ud.add_row(trans(X1.x));

//        if (norm(X0.x - X1.x) < epsilon || solution::f_calls > Nmax ||
//        solution::g_calls > Nmax) {
//            Xopt = X1;
//            X1.fit_fun(ff,ud1, ud2);
//            return X1;
//        }
//        X0 = X1;
//    }
            if (norm(X0.x - X1.x) < epsilon) {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if (std::max(solution::f_calls, solution::g_calls)) {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }
            X0 = X1;
        }
        return Xopt;
    }catch(string ex_info){
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution
CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0,
   double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        int n = get_len(x0);
        solution Xopt;
        solution X0, X1;
        X0.x = x0;
        matrix d(n, 1), P(n, 2);
        solution h;
        double *ab {};
        double beta;
//        Xopt.grad(gf);
//        d = -Xopt.g;
        d = -X0.grad(gf);
        while (true) {
            if (h0 < 0) {
//                P[0] = Xopt.x;
//                P[1] = d;
                P.set_col(X0.x, 0);
                P.set_col(d, 1);
                ab = expansion(ff,0, 1, 1.2, Nmax, ud1, P);
                h = golden(ff,ab[0], ab[1], epsilon, Nmax, ud1, P);
                X1.x = X0.x + h.x * d;
            } else {
                X1.x = X0.x + h0 * d;
            }
//            ud->add_row(trans(X1.x));
            Xopt.ud.add_row(trans(X1.x));
//            if (norm(Xopt1.x - Xopt.x) < epsilon || solution::f_calls > Nmax ||
//                solution::g_calls > Nmax) {
//                Xopt1.fit_fun(ff,ud1,ud2);
//                return Xopt1;
//            }
            if (norm(X1.x - X0.x) < epsilon) {
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
            if(std::max(solution::f_calls, solution::g_calls) > Nmax){
                Xopt = X1;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }

            X1.grad(gf);
            beta = pow(norm(X1.g), 2) / pow(norm(X0.g), 2);
            d = -X1.g + beta * d;
            X0 = X1;

//            Xopt1.grad(gf);
//             betda do okreslania kolejnych wartosci kierunku
//            beta = pow(norm(Xopt1.g), 2) / pow(norm(Xopt.g), 2);
//            d = -Xopt1.g + beta * d;
//            Xopt = Xopt1;
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
                matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        int n = get_len(x0);
        solution X0, X1;
        X0.x = x0;
        matrix d(n, 1);
        matrix P(n, 2);
        solution h;
        double *ab{};
        while (true) {
            X0.grad(gf);
            X0.hess(Hf);
            d = -inv(X0.H) * X0.g;
            if (h0 < 0) {
                P.set_col(X0.x, 0);
                P.set_col(d, 1);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, P);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, P);
                X1.x = X0.x + h.x * d;
            } else {
                X1.x = X0.x + h0 * d;
            }
            Xopt.ud.add_row(trans(X1.x));
            if (norm(X0.x - X1.x) < epsilon) {
                X1.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                return X1;
            }
            if (norm(X0.x - X1.x) < epsilon) {
                X1.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                return X1;
            }
            X0 = X1;
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution
golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1,
       matrix ud2) {
    try {
        solution Xopt;
        double alpha = (sqrt(5) - 1) / 2;
        solution A, B, C, D;
        A.x = a;
        B.x = b;
        C.x = B.x - alpha * (B.x - A.x);
        C.fit_fun(ff, ud1, ud2);
        D.x = A.x + alpha * (B.x - A.x);
        D.fit_fun(ff, ud1, ud2);
        while (true) {
            if (C.y < D.y) {
                B = D;
                D = C;
                C.x = B.x - alpha * (B.x - A.x);
                C.fit_fun(ff, ud1, ud2);
            } else {
                A = C;
                C = D;
                D.x = A.x + alpha * (B.x - A.x);
                D.fit_fun(ff, ud1, ud2);
            }
            if (B.x - A.x < epsilon) {
                Xopt.x = (A.x + B.x) / 2;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
//            A.x = (A.x + B.x) / 2;
//            A.fit_fun(ff, ud1, ud2);
//            return A;
            }
            if (std::max(solution::f_calls, solution::g_calls) > Nmax) {
                Xopt.x = (A.x + B.x) / 2;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
//            A.x = (A.x + B.x) / 2;
//            A.fit_fun(ff, ud1, ud2);
//            return A;
            }
        }
        return Xopt;
    }
    catch (string ex_info)
    {
        throw ("solution golden(...):\n" + ex_info);
    }
}

// #####################################################################################################################

solution
Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        int n = get_len(x0);
        matrix D = ident_mat(n);
        matrix A(n, 2);
        solution X, P, h;
        X.x = x0;
        double *ab{};
        while(true){
            P = X;
            for (int i = 0; i < n){
                A.set_col(P.x, 0);
                A.set_col(D[i], 1);
                ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
                h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
                P.x = P.x + h.x * D[i];
            }
            if(norm(X.x - P.x) < epsilon){
                Xopt = P;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 0;
                break;
            }
//            if(std::max(solution::f_calls, solution::g_calls) > Nmax){
            if(solution::f_calls > Nmax){
                Xopt = P;
                Xopt.fit_fun(ff, ud1, ud2);
                Xopt.flag = 1;
                break;
            }
            for (int i = 0; i < n - 1; i++){D.set_col(D[i+1], i);}
            D.set_col(P.x - X.x, n - 1);
            A.set_col(P.x, 0);
            A.set_col(D[n - 1], 1);
            ab = expansion(ff, 0, 1, 1.2, Nmax, ud1, A);
            h = golden(ff, ab[0], ab[1], epsilon, Nmax, ud1, A);
            X.x = P.x + h.x * D[n - 1];
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution
EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix limits, int mi, int lambda, matrix sigma0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        solution *P = new solution[mi + lambda];
        solution *Pm = new solution[mi];

        for(int i = 0; i < mi; i++){
            P[i].x = matrix(N, 2);
            for(int j = 0; j < N; j++){
                P[i].x(j, 0) = (limits(j, 1) - limits(j, 0)) * m2d(rand_mat(1, 1)) + limits(j, 0);
                P[i].x(j, 1) = sigma0(j);
            }
            if(P[i].fit_fun(ud1, ud2) < epsilon){
                Xopt = P[i];
                return Xopt;
            }
        }

        matrix IFF(mi, 1), t(N, 2);
        double r, S, s_IFF;

        while (true) {
            s_IFF = 0;
            for (int i = 0; i < mi; ++i) {
                IFF(i) = 1 / P[i].y();
                s_IFF += IFF(i);
            }
            for (int i = 0; i < lambda; ++i) {
                r = s_IFF * rand_mat(1, 1)();
                S = 0;
                for (int j = 0; j < mi; ++j) {
                    S += IFF(j);
                    if (r <= S) {
                        P[mi + i] = P[j];
                        break;
                    }
                }
            }

            for (int i = 0; i < lambda; ++i) {
                r = m2d(rand_mat(1, 1));
                for (int j = 0; j < N; ++j) {
                    P[mi + i].x(j, 1) *= exp(alfa * r + beta * m2d(rand_mat(1, 1));
                    P[mi + i].x(j, 0) += P[mi + i].x(j, 1) * m2d(rand_mat(1, 1));
                }
            }
            for (int i = 0; i < lambda; i += 2) {
                r = m2d(rand_mat(1, 1));
                t = P[mi + i].x;
                P[mi + i].x = r * P[mi + i].x + (1 - r) * P[mi + i + 1].x;
                P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r) * t;
            }
            if (solution::f_calls > Nmax) {
                return P[0];
            }
        }
        return Xopt;
    }

    catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}
