#include "payoff.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

double Asian_call_payoff(const std::vector<double>& path,double k)
{
    //const double sum = std::accumulate(path.begin()+1,path.end(),0.0);
    //const double average = sum / (path.size() - 1);
    if(path.size() < 2)
        throw std::invalid_argument("Path must contain at least two prices");
    double sum = 0.0;
    for(std::size_t i=1;i<path.size();i++)//debug here i should start with 1 rather than 0
        sum += path[i];
    const double average = sum / (path.size() - 1);
    return std::max(average - k,0.0);
}

double Asian_put_payoff(const std::vector<double>& path,double k)
{
    if(path.size() < 2)
        throw std::invalid_argument("Path must contain at least two prices");
    double sum = 0.0;
    for(std::size_t i=1;i<path.size();i++)
        sum += path[i];
    const double average = sum / (path.size() - 1);
    return std::max(k - average,0.0);
}
//Fixed-Strike Lookback Option
double Lookback_call_payoff(const std::vector<double>& path,double k)
{
    //auto it = std::max_element(path.begin(),path.end());
    //double max_price = *it;//Here, * is a dereference, not a pointer.
    //Returns the position (iterator) of the largest element rather than the element itself
    if (path.empty())
        throw std::invalid_argument("Path must not be empty");
    double max_price = path[0];
    for(std::size_t i=0;i<path.size();i++)
    {
        if(path[i] > max_price)
            max_price = path[i];
    }
    return std::max(max_price - k,0.0);
}
double Lookback_put_payoff(const std::vector<double>& path,double k)
{
    if (path.empty())
        throw std::invalid_argument("Path must not be empty");
    double min_price = path[0];
    for(std::size_t i=0;i<path.size();i++)
    {
        if(path[i] < min_price)
            min_price = path[i];
    }
    return std::max(k - min_price,0.0);
}
