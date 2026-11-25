#include "gbm.hpp"
#include "gbm_validation.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

std::vector<double> simulation_GBM(double s0,double sigma,double r,std::size_t steps,double T,std::mt19937_64 & rng)//rng =: Random Number Generator
//Exact Discretization not Euler–Maruyama!
{
    validate_gbm(s0, sigma, r, steps, T);
    const double dt = T / steps;//deta t
    std::normal_distribution<double> normal(0.0,1.0);
    std::vector<double> S(steps+1);//t0----tn
    S[0] = s0;
    for(size_t i=0;i<steps;i++)//debug i<step rather than r<S.size()
    {
        const double Z = normal(rng);//Generate a random number Z that follows a normal distribution.
        const double dirft = (r - 0.5*sigma*sigma)*dt;
        const double diffusion = sigma * std::sqrt(dt) * Z;
        S[i+1] = S[i] * std::exp(dirft + diffusion);
    }
    return S;// return an array
}
