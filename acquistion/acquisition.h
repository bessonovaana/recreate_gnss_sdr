
#ifndef ACQUISITION_H
#define ACQUISITION_H

#include <complex>
#include <cstdint>
#include <vector>

struct AcquisitionResult
{
    float doppler = 0.0f;
    long long delay = 0;
    float glrt = 0.0f;
    bool detected = false;
};

// Генерация локального кода ГЛОНАСС
std::vector<std::complex<float>> generate_glonass_code(
    int32_t fs,
    uint32_t chip_shift
);

// Захват сигнала
AcquisitionResult acquistion(
    const std::vector<std::complex<float>>& iq_samples,
    const std::vector<std::complex<float>>& local_code,
    float fs,
    float f_min,
    float f_max,
    float f_step,
    float threshold
);

#endif
