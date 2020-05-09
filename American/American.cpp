#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <algorithm>
#include <cassert>

// S: spot price at time 0
// T: time to maturity (in years)
// v: the volatility
// r: risk-free rate
// q: dividend
// n: the number of time steps
// K: the strike of the american call option

enum OptionType { Call, Put };

//Pricing American option with binomial method
double amer(OptionType optType, double K, double T, double S, double sigma, double r, double q, int n)
{
    std::vector<double> states(n + 1);
    std::vector<double> payoff(n + 1);
    double dt = T / n;
    double u = exp(sigma * sqrt(dt));
    double p = (std::exp((r - q) * dt) - 1 / u) / (u - 1 / u);
    // initialize the final states, apply payoff directly
    for (int i = 0; i <= n; i++) {
        double Si = S * std::pow(u, n - 2 * i);
        switch (optType) {
        case Call:
            states[i] = Si > K ? Si - K : 0;
        case Put:
            states[i] = Si < K ? K - Si : 0;       
        }
        states[i] = Si < K ? K - Si : 0;
    }
    for (int k = n - 1; k >= 0; k--) {
        for (int i = 0; i <= k; i++) {
            switch (optType)
            {
            case Call:
                payoff[i] = S * pow(u, k - 2 * i) - K;
                break;
            case Put:
                payoff[i] = K - S * pow(u, k - 2 * i);
                break;
            default:
                throw "unsupported optionType";
            } 
            states[i] = exp(-r * dt) * (states[i] * p + states[i + 1] * (1 - p));
            states[i] = std::max(states[i], payoff[i]);
        }
    }
    return states[0];
}

//pricing American option value with the smooth Black-Scholes analytic formula
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

double amerSmooth(OptionType optType,double S, double K, double T, double sigma, double r, double q, double n)
{

    std::vector<double> states(n);
    std::vector<double> euroValue(n);
    std::vector<double> payoff(n);
    double dt = T / n;
    double u = exp(sigma * sqrt(dt));
    double p = (std::exp((r-q) * dt) - 1 / u) / (u - 1 / u);

    for (int i = 0; i <= n - 1; i++) {
        double Si = S * std::pow(u, n - 1 - 2 * i);
        // initialize the second last payoff directly
        switch (optType)
        {
        case Call:
            payoff[i] = Si > K ? Si - K : 0;
            break;
        case Put:
            payoff[i] = Si < K ? K - Si : 0;
            break;
        default:
            throw "unsupported optionType";
        }
        // initialize the second last continuation
        euroValue[i] = BlackScholesEuro(optType, K, dt, Si, sigma, r, q);
        states[i] = std::max(payoff[i], euroValue[i]);
    }
    for (int k = n - 2; k >= 0; k--) {
        for (int i = 0; i <= k; i++) {
            switch (optType)
            {
            case Call:
                payoff[i] = S * pow(u, k - 2 * i) - K;
                break;
            case Put:
                payoff[i] = K - S * pow(u, k - 2 * i);
                break;
            default:
                throw "unsupported optionType";
            }            
            states[i] = exp(-r * dt) * (states[i] * p + states[i + 1] * (1 - p));
            states[i] = std::max(states[i], payoff[i]);
        }
    }
    return states[0];
}

//Under Black - Scholes model, one way to approximate American put option price is the
//quadratic approximation method by Barone - Adesi and Whaley(1987).
bool RootBracketing(std::function<double(double)> f, double& a, double& b
    , double min = std::numeric_limits<double>::min()
    , double max = std::numeric_limits<double>::max())
{
    const int NTRY = 50;
    const double FACTOR = 1.6;
    if (a >= b) throw("wrong input a and b in RootBracketing");
    double f1 = f(a);
    double f2 = f(b);
    for (int j = 0; j < NTRY; j++) {
        if (f1 * f2 < 0.0) return true;
        if (std::abs(f1) < std::abs(f2))
            f1 = f(a = std::max(a + FACTOR * (a - b), min));
        else
            f2 = f(b = std::min(b + FACTOR * (b - a), max));
    }
    return false;
}

double rfbisect(std::function<double(double)> f, double& a, double& b, double tol)
{
    assert(a < b && f(a) * f(b) < 0);
    double c;
    while ((b - a) / 2 > tol) {
        c = (a + b) / 2.0;
        //std::cout << "(a, b) = (" << a << ", " << b << ")" << std::endl;
        if (std::abs(f(c)) < tol)
            return c;
        else {
            if (f(a) * f(c) < 0)
                b = c;
            else
                a = c;
        }
    }
    return c;
}


double whaley(OptionType optType, double S_0, double K, double T, double sigma, double r, double q)
{
    double b = r - q;
    double alpha = 2 * r / pow(sigma, 2);
    double beta = 2 * b / pow(sigma, 2);
    //root searcher to get S_hat
    double upperbound = S_0, lowerbound = 0;
    auto fun = [optType,T, r, sigma, K,q,b](double S) {
        double M = 1 - 2 * b / pow(sigma, 2);
        double d1_hat = (std::log(S / K) + (b + 0.5 * pow(sigma, 2)) * T) / (sigma * sqrt(T));
        double q_1 = (M - std::sqrt(pow(M, 2) + 8 * r / pow(sigma, 2) / (1 - std::exp(-r * T)))) / 2;
        double q_2= (M + std::sqrt(pow(M, 2) + 8 * r / pow(sigma, 2) / (1 - std::exp(-r * T)))) / 2;
        double A_1 = -(S / q_1) * (1 - exp(-q*T)*cnorm(-d1_hat));
        double A_2= (S / q_2) * (1 - exp(-q * T) * cnorm(d1_hat));
        double V_hat;
        double V;
        switch (optType)
        {
        case Call:
            V_hat = BlackScholesEuro(Call, K, T, S, sigma, r, q);
            V = V_hat + A_2 - S + K;
            break;
        case Put:
            V_hat = BlackScholesEuro(Put, K, T, S, sigma, r, q);
            V = V_hat + A_1 + S - K;
            break;
        default:
            throw "unsupported optionType";
        }  
        return V;
    };
    RootBracketing(fun, lowerbound, upperbound,1e-6);
    double S_hat = rfbisect(fun, lowerbound, upperbound, 1e-6);
    //get option value
    double M = 1 - 2 * b / pow(sigma, 2);
    double d1_hat = (std::log(S_hat / K) + (b + 0.5 * pow(sigma, 2)) * T) / (sigma * sqrt(T));
    double q_1 = (M - std::sqrt(pow(-M, 2) + 8 * r / (pow(sigma, 2) * (1 - std::exp(-r * T))))) / 2;
    double q_2 = (M + std::sqrt(pow(-M, 2) + 8 * r / (pow(sigma, 2) * (1 - std::exp(-r * T))))) / 2;
    double A_1 = -(S_hat / q_1) * (1 - exp(-q * T) * cnorm(-d1_hat));
    double A_2 = (S_hat / q_2) * (1 - exp(-q * T) * cnorm(d1_hat));
    double V_0;
    switch (optType)
    {
    case Call:
        if (S_0 < S_hat) {
            V_0 = BlackScholesEuro(Call, K, T, S_0, sigma, r, q) + A_2 * pow((S_0 / S_hat), q_2);
        }
        else {
            V_0 = S_0 - K;
        }
        break;
    case Put:
        if (S_0 > S_hat) {
            V_0 = BlackScholesEuro(Put, K, T, S_0, sigma, r, q) + A_1 * pow((S_0 / S_hat), q_1);
        }
        else {
            V_0 = K - S_0;
        }
        break;
    default:
        throw "unsupported optionType";
    }
    
    return  V_0;
}


int main() {
    double S, K, T, sigma, r, tenors, q;
    S = 100.0; K = 90.0; T = 1.0; sigma = 0.15; r = 0.03, q = 0.01, tenors = 100;
    std::cout << amer(Put, K, T, S, sigma, r, q, tenors) << std::endl;
    std::cout << amerSmooth(Put, S, K, T, sigma, r, q, tenors) << std::endl;
    std::cout << whaley(Put,S, K, T, sigma, r, q) << std::endl;
}
