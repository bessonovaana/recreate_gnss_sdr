
#ifndef IQ_SAMPLES_H
#define IQ_SAMPLES_H

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

struct IQSamples
{
    std::vector<std::complex<float>> samples;
    double sampling_frequency_hz = 0.0;

    bool load(const std::string& filename);
    bool save(const std::string& filename) const;

    std::size_t size() const;
    bool empty() const;
};

#endif
