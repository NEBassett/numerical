#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <iostream>


auto gramPolynomial(const unsigned deg, const double x, const int N)
{
  const auto factorial = [](int x){
    double prod = 1.0f;
    for(;x > 0;x--)
    {
      prod = prod*x;
    }
    return prod;
  };

  const auto factorialPower = [](double x, double n){
    return std::tgamma(x+1)/std::tgamma(x-n+1);
  };

  double sum = 0.0f;
  for(size_t i = 0; i <= deg; i++)
  {
    sum = sum + pow(-1,i)*(factorialPower(deg+i,2*i)/pow(factorial(i),2))*(factorialPower(x,i)/factorialPower(N,i));
  }
  return sum*pow(-1,deg);
}

template<int r>
class gramApproximationPolynomial
{
  const int N;
  std::vector<double> targets;
  std::vector<double> coefficients;
public:
  auto operator()(double x)
  {
    double sum = 0.0f;

    for(size_t i = 0; i <= r; i++)
    {
      sum += coefficients[i]*gramPolynomial(i, x, N);
    }

    return sum;
  }

  gramApproximationPolynomial(std::initializer_list<double> values, double a, double b) : N(values.size() - 1), targets(values)
  {
    const auto factorial = [](int x){
      double prod = 1.0f;
      for(;x > 0;x--)
      {
        prod = prod*x;
      }
      return prod;
    };

    coefficients.reserve(r+1);

    if(N%2 != 0 || a>=b)
    {
      throw std::domain_error("Number of ordinates must be odd and range must be valid.\n");
    }

    const double h = (b-a)*(1/N);

    for(size_t i = 0; i <= r; i++)
    {
      double yR = 0.0f;
      double aR = 0.0f;

      for(size_t j = -N/2; j != (N/2 + 1); j++)
      {
        auto val = gramPolynomial(i, N/2 + j, N);
        yR += pow(val, 2);
        aR += val*targets[j+N/2];
      }

      coefficients[i] = aR/yR;
    }
  }
};

int main()
{
  gramApproximationPolynomial<4> poly({1.0f,4.0f,2.60f,4.0f,5.67f}, 1.0f, 2.0f);
  std::cout << poly(1.745f) << std::endl;
}
