#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <vector>

const uint32_t L =511;
const int fc = L * 1000;

std::vector<std::complex<float>>  generate_glonass_code (int32_t fs, uint32_t chip_shift);

std::vector<float> correlation_power(const std::vector<std::complex<float>>& correlation);




