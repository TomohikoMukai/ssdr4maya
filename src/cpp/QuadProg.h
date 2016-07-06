#ifndef QUADPROG_H
#define QUADPROG_H
#pragma once

#include <Eigen/Core>
#include "QuadProg++.hh"

// min (0.5 * xv^T * gm * xv + gv^T * x)
//  s.t. cem * xv + cev = 0
//       cim * xv + civ >= 0
double SolveQP(const Eigen::MatrixXd& gm, const Eigen::VectorXd& gv,
    const Eigen::MatrixXd& cem, const Eigen::VectorXd& cev,
    const Eigen::MatrixXd& cim, const Eigen::VectorXd& civ,
    Eigen::VectorXd& xv);
double TestSolveQP();

#endif //QUADPROG_H
