#include "QuadProg.h"

using namespace Eigen;

double SolveQP(const MatrixXd& gm, const VectorXd& gv,
    const MatrixXd& cem, const VectorXd& cev,
    const MatrixXd& cim, const VectorXd& civ,
    VectorXd& xv)
{
    QuadProgPP::Matrix<double> G(gm.rows(), gm.cols());
    for (int i = 0; i < gm.rows(); ++i)
    {
        for (int j = 0; j < gm.cols(); ++j)
        {
            G[i][j] = gm(i, j);
        }
    }
    QuadProgPP::Vector<double> g0(gv.size());
    for (int i = 0; i < gv.size(); ++i)
    {
        g0[i] = gv[i];
    }
    QuadProgPP::Matrix<double> CE;
    if (cem.rows() > 0 && cem.cols() > 0)
    {
        CE.resize(cem.cols(), cem.rows());
        for (int i = 0; i < cem.rows(); ++i)
        {
            for (int j = 0; j < cem.cols(); ++j)
            {
                CE[j][i] = cem(i, j);
            }
        }
    }
    QuadProgPP::Vector<double> ce0;
    if (cev.size() > 0)
    {
        ce0.resize(cev.size());
        for (int i = 0; i < cev.size(); ++i)
        {
            ce0[i] = cev[i];
        }
    }
    QuadProgPP::Matrix<double> CI;
    if (cim.rows() > 0 && cim.cols() > 0)
    {
        CI.resize(cim.cols(), cim.rows());
        for (int i = 0; i < cim.rows(); ++i)
        {
            for (int j = 0; j < cim.cols(); ++j)
            {
                CI[j][i] = cim(i, j);
            }
        }
    }
    QuadProgPP::Vector<double> ci0;
    if (civ.size() > 0)
    {
        ci0.resize(civ.size());
        for (int i = 0; i < civ.size(); ++i)
        {
            ci0[i] = civ[i];
        }
    }
    QuadProgPP::Vector<double> x(0.0, xv.size());
    double retval = QuadProgPP::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    for (int i = 0; i < xv.size(); ++i)
    {
        xv[i] = x[i];
    }
    return retval;
}

double TestSolveQP()
{
    MatrixXd H(2, 2);
    H(0, 0) = 1.0;
    H(0, 1) = -1.0;
    H(1, 0) = -1.0;
    H(1, 1) = 2.0;

    VectorXd f(2);
    f(0) = -2.0;
    f(1) = -6.0;

    MatrixXd A(2, 5);
    A(0, 0) = -1.0;
    A(1, 0) = -1.0;
    A(0, 1) = 1.0;
    A(1, 1) = -2.0;
    A(0, 2) = -2.0;
    A(1, 2) = -1.0;
    A(0, 3) = 1.0;
    A(1, 3) = 0.0;
    A(0, 4) = 0.0;
    A(1, 4) = 1.0;

    VectorXd b(5);
    b(0) = 2.0;
    b(1) = 2.0;
    b(2) = 3.0;
    b(3) = 0;
    b(4) = 0;

    MatrixXd C;
    VectorXd d;

    VectorXd xv(2);
    SolveQP(H, f, C, d, A, b, xv);
    printf("%f - %f\n", xv[0], xv[1]);
    return (xv[0] - 2.0 / 3.0) * (xv[0] - 2.0 / 3.0) + (xv[1] - 4.0 / 3.0) * (xv[1] - 4.0 / 3.0);
}
