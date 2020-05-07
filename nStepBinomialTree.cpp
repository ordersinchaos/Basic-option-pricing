#include <iostream>
#include <vector>
#include <cmath> // for std::exp()
#include <cassert> // for assertion on inputs
enum OptionType {Call, Put};

double nStepBinomialTree ( OptionType optType, double K, double T
			 , double S_0, double sigma, double rate, int N)
{
  std::vector<double> states(N+1);
  double dt = T / N;
  double b = std::exp((2*rate+sigma*sigma)*dt)+1;
  double u = (b + std::sqrt(b*b - 4*std::exp(2*rate*dt))) / 2 / std::exp(rate*dt);
  double p = (std::exp(rate*dt) -  1/u) / (u - 1/u);
  // initialize the final states, apply payoff directly
  for (int i = 0; i <= N; i++) {
    double S = S_0 * std::pow(u, N-2*i);
    switch (optType)
      {
      case Call:
	states[i] = S > K ? S - K : 0;
	break;
      case Put:
	states[i] = S < K ? K - S : 0;
	break;
      default:
	throw "unsupported optionType";
      }
  }
  for (int k = N-1; k >= 0; k--) 
    for (int i = 0; i <= k; i++)
      states[i] = states[i]*p + states[i+1]*(1-p);
  
  return exp(-rate*T) * states[0];
}
double cnorm(double);
double bsPricer(OptionType optType, double K, double T, double S_0, double sigma, double rate)
{
  double sigmaSqrtT = sigma * std::sqrt(T);
  double d1 = (std::log(S_0 / K) + rate)/sigmaSqrtT + 0.5 * sigmaSqrtT;
  double d2 = d1 - sigmaSqrtT;

  double V_0;
  switch (optType)
    {
    case Call:
      V_0 = S_0 * cnorm(d1) - K * exp(-rate*T) * cnorm(d2);
      break;
    case Put:
      V_0 = K * exp(-rate*T) * cnorm(-d2) - S_0 * cnorm(-d1);
      break;
    default:
      throw "unsupported optionType";
    }
  return V_0;
}

double cnorm(double x)
{
  // constants
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x)/sqrt(2.0);
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  return 0.5*(1.0 + sign*y);
}     

int main()
{
  OptionType optType;
  double K, T, S_0, S_u, S_d, rate, V_u, V_d, sigma;
  std::cout << "option type (0 for Call and 1 for Put): ";
  int optType_;
  std::cin >> optType_;
  assert(optType_ == 0 || optType_ == 1);
  optType = static_cast<OptionType>(optType_);
  std::cout << "strke: ";
  std::cin >> K;
  std::cout << "time to maturity: ";
  std::cin >> T;
  std::cout << "current value of the underlying stock: ";
  std::cin >> S_0;
  std::cout << "risk free interest rate: ";
  std::cin >> rate;
  std::cout << "volatility: ";
  std::cin >> sigma;

  double V_bs = bsPricer(optType, K, T, S_0, sigma, rate);
  std::cout << "The option price of black-scholes formula is "
	    << std::endl << V_bs << std::endl;
  for (int i = 1; i < 200; i++)
    std::cout << i << " " << nStepBinomialTree(optType, K, T, S_0, sigma, rate, i) - V_bs << std::endl;

  return 0;
}
