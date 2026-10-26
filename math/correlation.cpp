#include <math_func.h>

std::vector<float> correlation_power(const std::vector<std::complex<float>>& correlation){
    std::vector<float> power(correlation.size());
    for (std::size_t k = 0; k < correlation.size(); k++){
        power[k] = std::norm(correlation[k]);
    }
    return power;
}