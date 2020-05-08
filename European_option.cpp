#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <utility>
#include <algorithm>
#include <random>

enum OptionType {Call,Put};
double cnorm(double x)
{
    // constants
    double a1 = 0.254829592;
    double a2 = -0.284496736;
    double a3 = 1.421413741;
    double a4 = -1.453152027;
    double a5 = 1.061405429;
    double p = 0.3275911;
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x) / sqrt(2.0);
    double t = 1.0 / (1.0 + p * x);
    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);
    return 0.5 * (1.0 + sign * y);
}

double BlackScholesEuro(OptionType optType, double K, double T, double S_0, double sigma, double rate, double q)
{
    double sigmaSqrtT = sigma * std::sqrt(T);
    double d1 = (std::log(S_0 / K) + (rate - q) * T) / sigmaSqrtT + 0.5 * sigmaSqrtT;
    double d2 = d1 - sigmaSqrtT;

    double V_0;
    switch (optType)
    {
    case Call:
        V_0 = S_0 * exp(-q * T) * cnorm(d1) - K * exp(-rate * T) * cnorm(d2);
        break;
    case Put:
        V_0 = K * exp(-rate * T) * cnorm(-d2) - S_0 * exp(-q * T) * cnorm(-d1);
        break;
    default:
        throw "unsupported optionType";
    }
    return V_0;
}

double BinomialTreeEuro(OptionType optType, double K, double T, double S_0, double sigma, double rate, double q, int tenors) 
{
    std::vector<double> states(tenors + 1);
    double dt = T / tenors;
    double u = std::exp(sigma * sqrt(dt));
    double p = (std::exp((rate-q) * dt) - 1 / u) / (u - 1 / u);
    // initialize the final states, apply payoff directly
    for (int i = 0; i <= tenors; i++) {
        double S = S_0 * std::pow(u, tenors - 2 * i);
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
    for (int k = tenors - 1; k >= 0; k--)
        for (int i = 0; i <= k; i++)
            states[i] = exp(-rate * dt) * (states[i] * p + states[i + 1] * (1 - p));

    return states[0];
}

std::pair<double, double> MonteCaeloEuro(OptionType optType, int nPaths
    , double K, double S, double T, double rate, double q, double sigma)
{
    // Generate a normal distribution around that mean
    std::seed_seq seed{ 5 };
    std::mt19937 e(seed);
    std::normal_distribution<> normal_dist(0, 1);
    double sum = 0;
    double hsquare = 0;
    for (int i = 0; i < nPaths; i++) {
        double wT = normal_dist(e) * sqrt(T);
        double h = S * std::exp((rate - q - 0.5 * sigma * sigma) * T + sigma * wT);
        switch (optType) {
        case Call:
            h = h > K ? h - K : 0;
            break;
        case Put:
            h = h < K ? K - h : 0;
            break;
        default:
            throw "unsupported optionType";
        }      
        sum += h;
        hsquare += h * h;
    }
    double pv = std::exp(-rate * T) * sum / nPaths;
    double Stderr = sqrt((hsquare / nPaths - (sum / nPaths) * (sum / nPaths)) / nPaths);
    return std::pair<double, double>(pv, Stderr);
}

int main() {
    double K, T, S_0, sigma, rate, q;
    int tenors,nPath;
    K = 90.0, T = 1.0, S_0 = 100.0, sigma = 0.15, rate = 0.03, q = 0.01, tenors = 1000, nPath = 100000;
    std::cout << BlackScholesEuro(Call, K, T, S_0, sigma, rate, q) << std::endl;
    std::cout << BinomialTreeEuro(Call, K, T, S_0, sigma, rate, q, tenors) << std::endl;
    std::cout << MonteCaeloEuro(Call, nPath, K, S_0, T, rate, q, sigma).first << std::endl;
}
