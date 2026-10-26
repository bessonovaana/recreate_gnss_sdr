#ifndef TRACKING_H
#define TRACKING_H

#include <array>
#include <complex>
#include <vector>

#include "acquisition.h"

struct TrackingState
{
    // Текущее состояние генератора кода.
    double code_phase_chips = 0.0;
    double code_rate_chips_s = 511000.0;

    // Текущее состояние генератора несущей.
    double carrier_phase_rad = 0.0;
    double carrier_doppler_hz = 0.0;

    // Ошибки дискриминаторов.
    double dll_error = 0.0;
    double pll_error = 0.0;

    // Номинальные параметры, полученные из Acquisition.
    double code_rate_nominal_chips_s = 511000.0;
    double carrier_doppler_nominal_hz = 0.0;

    // Состояния интеграторов DLL и PLL.
    double dll_integrator_chips_s = 0.0;
    double pll_integrator_rad_s = 0.0;

    // Состояние слежения.
    bool loops_initialized = false;
    bool locked = false;
};

struct TrackingCorrelatorResult
{
    std::complex<float> early{0.0f, 0.0f};
    std::complex<float> prompt{0.0f, 0.0f};
    std::complex<float> late{0.0f, 0.0f};
};

TrackingCorrelatorResult tracking_correlate(
    const std::vector<std::complex<float>>& iq_samples,
    const std::array<int, 511>& code,
    TrackingState& state,
    double fs
);

double dll(const TrackingCorrelatorResult& corr);

double pll(const TrackingCorrelatorResult& corr);

TrackingCorrelatorResult tracking_process_block(
    const std::vector<std::complex<float>>& iq_samples,
    const std::array<int, 511>& code,
    TrackingState& state,
    double fs
);

#endif
