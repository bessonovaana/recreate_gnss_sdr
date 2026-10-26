#include "iq_samples.h"

#include <fstream>
#include <iostream>


bool IQSamples::load(const std::string& filename)
{
    std::ifstream file(filename, std::ios::binary);

    if (!file)
    {
        std::cerr << "Cannot open file: " << filename << '\n';
        return false;
    }

    samples.clear();

    float i_value;
    float q_value;

    while (true)
    {
        if (!file.read(
                reinterpret_cast<char*>(&i_value), sizeof(float)))
        {
            break;
        }

        if (!file.read(
                reinterpret_cast<char*>(&q_value), sizeof(float)))
        {
            std::cerr
                << "Warning: incomplete IQ sample at end of file. "
                << "Ignoring the final incomplete sample.\n";
            break;
        }

        samples.emplace_back(i_value, q_value);
    }

    if (samples.empty())
    {
        std::cerr << "No complete IQ samples found\n";
        return false;
    }

    if (file.bad())
    {
        std::cerr << "Error while reading file\n";
        samples.clear();
        return false;
    }

    return true;
}

bool IQSamples::save(const std::string& filename) const
{
    std::ofstream file(
        filename, std::ios::binary | std::ios::trunc);

    if (!file)
    {
        std::cerr << "Cannot create file: " << filename << '\n';
        return false;
    }

    for (const auto& sample : samples)
    {
        float i_value = sample.real();
        float q_value = sample.imag();

        file.write(
            reinterpret_cast<const char*>(&i_value), sizeof(float));

        file.write(
            reinterpret_cast<const char*>(&q_value), sizeof(float));

        if (!file)
        {
            std::cerr << "Error while writing IQ samples\n";
            return false;
        }
    }

    return true;
}

std::size_t IQSamples::size() const
{
    return samples.size();
}

bool IQSamples::empty() const
{
    return samples.empty();
}
