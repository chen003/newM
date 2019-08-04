#include "eigen_calc.h"

#include <algorithm>
#include <ccomplex>
#include <vector>

#include <eigen3/Eigen/Eigenvalues>

//#include <iostream>

extern "C" void GetEigenValuesAndVectors(int n, double *matData,
                                         double _Complex *eigenValueOut,
                                         double _Complex *eigenVectorOut) {
  struct {
    double *data;
    Eigen::Index rows, cols;
  } cMatrix = {matData, n, n};
  const Eigen::MatrixXd &eigenMatrix =
      reinterpret_cast<Eigen::MatrixXd &>(cMatrix);
  Eigen::EigenSolver<Eigen::MatrixXd> es(eigenMatrix, true);

  auto eigenValues = es.eigenvalues();
  auto eigenVectors = es.eigenvectors();

  //    for (int i = 0; i < n; i++) {
  //        std::cout << eigenValues(i, 0) << '\t';
  //        for (int j = 0; j < n; j++) {
  //            std::cout << eigenVectors(j, i) << ' ';
  //        }
  //        std::cout << "\n\n";
  //    }

  struct EigenValueAndVector {
    std::complex<double> EigenValue;
    std::vector<std::complex<double>> EigenVector;

    EigenValueAndVector(std::complex<double> eValue,
                        decltype(eigenVectors.col(0)) &&eVector) {
      EigenValue = eValue;
      for (int i = 0; i < eVector.rows(); i++)
        EigenVector.push_back(eVector(i, 0));
    }
  };

  std::vector<EigenValueAndVector> eigenValuesAndVectors;
  for (int i = 0; i < n; i++)
    eigenValuesAndVectors.emplace_back(eigenValues(i, 0), eigenVectors.col(i));

  auto realPartLess = [](const EigenValueAndVector &x,
                         const EigenValueAndVector &y) {
    return x.EigenValue.real() < y.EigenValue.real();
  };
  std::sort(eigenValuesAndVectors.begin(), eigenValuesAndVectors.end(),
            realPartLess);

  for (int i = 0; i < n; i++) {
    eigenValueOut[i] = reinterpret_cast<const double _Complex &>(
        eigenValuesAndVectors[i].EigenValue);
    std::copy(eigenValuesAndVectors[i].EigenVector.begin(),
              eigenValuesAndVectors[i].EigenVector.end(),
              reinterpret_cast<std::complex<double> *>(eigenVectorOut + n * i));
  }
}
