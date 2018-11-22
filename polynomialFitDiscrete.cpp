#include <Eigen/Dense>
#include <cmath>
#include <iostream>

template<int N>
class fittedPolynomial{
  Eigen::Matrix<float, 1, N> coefficients;
public:
  auto operator()(float x)
  {
    float value = 0;
    for(size_t i = 0; i < N; i++)
    {
      value = value + coefficients(0,i)*pow(x, i);
    }
    return value;
  }


  fittedPolynomial(Eigen::Matrix<float, 2, N> &values)
  {
    Eigen::Matrix<float, N, N> sys;
    Eigen::Matrix<float, 1, N> results;
    results(0) = values(1,0);
    results(1) = values(1,1);
    results(2) = values(1,2);

    for(size_t i = 0; i < N; i++)
    {
      for(size_t j = 0; j < N; j++)
      {
        sys(i,j) = pow(values(0, i), j);
      }
    }

    coefficients = sys.colPivHouseholderQr().solve(results.transpose());
  }
};

int main()
{
  Eigen::Matrix<float, 2, 3> targets(2,3);
  targets(0,0) = 3;
  targets(0,1) = 8;
  targets(0,2) = 23;
  targets(1,0) = 24;
  targets(1,1) = 5;
  targets(1,2) = -3;

  fittedPolynomial<3> poly(targets);

  std::cout << poly(3) << ", " << poly(8)<< ", " << poly(23);
}
