#include "tracking.h"
#include "math_func.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

namespace
{
    constexpr double PI = 3.14159265358979323846;
    constexpr double TWO_PI = 2.0 * PI;

    // Смещение Early и Late относительно Prompt, в чипах.
    constexpr double EARLY_LATE_SPACING = 0.5;

    // Параметры петлевых фильтров.
    // Это начальные значения для экспериментов, а не универсальные настройки.
    constexpr double DLL_BANDWIDTH_HZ = 1.0;
    constexpr double PLL_BANDWIDTH_HZ = 15.0;
    constexpr double LOOP_DAMPING = 0.7071067811865476;

    int wrap_index(int index)
    {
        index %= L;

        if (index < 0)
        {
            index += L;
        }

        return index;
    }

    int code_index(double phase)
    {
        return wrap_index(
            static_cast<int>(std::floor(phase)));
    }

    double wrap_code_phase(double phase)
    {
        phase = std::fmod(phase, static_cast<double>(L));

        if (phase < 0.0)
        {
            phase += L;
        }

        return phase;
    }

    // Натуральная частота петлевого фильтра.
    double loop_wn(double bandwidth_hz)
    {
        return bandwidth_hz * 8.0 * LOOP_DAMPING / (4.0 * LOOP_DAMPING * LOOP_DAMPING + 1.0);
    }

    // Фазовый дискриминатор Costas для BPSK.
    // Устойчив к смене знака навигационного символа.
    double costas_pll(const TrackingCorrelatorResult& corr)
    {
        const double i = static_cast<double>(corr.prompt.real());
        const double q = static_cast<double>(corr.prompt.imag());

        if (i == 0.0 && q == 0.0)
        {
            return 0.0;
        }

        return std::atan2(
            (i >= 0.0) ? q : -q,
            std::abs(i));
    }
}


 
// Корреляция Early / Prompt / Late
 

TrackingCorrelatorResult tracking_correlate(const std::vector<std::complex<float>>& iq_samples,const std::array<int, L>& code,TrackingState& state,double fs)
{
    TrackingCorrelatorResult result{};

    if (fs <= 0.0 || !std::isfinite(fs) || iq_samples.empty())
    {
        return result;
    }

    if (!std::isfinite(state.code_phase_chips) ||
        !std::isfinite(state.code_rate_chips_s) ||
        !std::isfinite(state.carrier_phase_rad) ||
        !std::isfinite(state.carrier_doppler_hz))
    {
        return result;
    }

    // Скорость изменения фазы кода на один отсчёт.
    const double code_step = state.code_rate_chips_s / fs;

    // Скорость изменения фазы несущей на один отсчёт.
    const double carrier_step =TWO_PI * state.carrier_doppler_hz / fs;

    for (const auto& sample : iq_samples)
    {
        const double prompt_phase = state.code_phase_chips;
        const int early_index = code_index(prompt_phase - EARLY_LATE_SPACING);
        const int prompt_index = code_index(prompt_phase);
        const int late_index = code_index( prompt_phase + EARLY_LATE_SPACING);

        // Комплексная локальная несущая с противоположной фазой.
        const std::complex<float> carrier(static_cast<float>(std::cos(state.carrier_phase_rad)),static_cast<float>(-std::sin(state.carrier_phase_rad)));

        // Компенсация несущей.
        const std::complex<float> wiped = sample * carrier;

        // Корреляция с тремя локальными копиями кода.
        result.early +=wiped * static_cast<float>(code[early_index]);
        result.prompt +=wiped * static_cast<float>(code[prompt_index]);
        result.late +=wiped * static_cast<float>(code[late_index]);

        // Сохраняем непрерывность фазы несущей.
        state.carrier_phase_rad = std::remainder(state.carrier_phase_rad + carrier_step,TWO_PI);

        // Сохраняем непрерывность фазы кода.
        state.code_phase_chips = wrap_code_phase(state.code_phase_chips + code_step);
    }

    return result;
}


 
// DLL-дискриминатор
 

double dll(const TrackingCorrelatorResult& corr)
{
    const double early = std::abs(corr.early);
    const double late = std::abs(corr.late);

    const double denominator = early + late;

    if (denominator <= 0.0 || !std::isfinite(denominator))
    {
        return 0.0;
    }

    return (early - late) / denominator;
}


 
// PLL-дискриминатор
 

double pll(const TrackingCorrelatorResult& corr)
{
    return costas_pll(corr);
}


 
// Обработка блока и обновление петель слежения
TrackingCorrelatorResult tracking_process_block(const std::vector<std::complex<float>>& iq_samples,const std::array<int, L>& code,TrackingState& state,double fs)
{
    TrackingCorrelatorResult corr{};

    if (fs <= 0.0 || !std::isfinite(fs) || iq_samples.empty())
    {
        return corr;
    }

    const double block_time =
        static_cast<double>(iq_samples.size()) / fs;

    // Первый вызов: запоминаем начальные параметры,
    // полученные, например, из Acquisition.
    if (!state.loops_initialized)
    {
        state.code_rate_nominal_chips_s =state.code_rate_chips_s;

        state.carrier_doppler_nominal_hz =state.carrier_doppler_hz;

        state.dll_integrator_chips_s = 0.0;
        state.pll_integrator_rad_s = 0.0;

        state.loops_initialized = true;
    }

    // 1. Корреляция входного сигнала с локальным кодом.
    corr = tracking_correlate(iq_samples, code, state, fs);

    // 2. Вычисление ошибок.
    state.dll_error = dll(corr);
    state.pll_error = pll(corr);

    if (!std::isfinite(state.dll_error) ||
        !std::isfinite(state.pll_error))
    {
        return corr;
    }

    // 3. Параметры DLL.
    const double dll_wn = loop_wn(DLL_BANDWIDTH_HZ);
    const double dll_k1 =2.0 * LOOP_DAMPING * dll_wn;
    const double dll_k2 =dll_wn * dll_wn;

    // DLL-дискриминатор нормирован.
    // Учитываем выбранное смещение Early/Late.
    const double dll_phase_scale = EARLY_LATE_SPACING;

    // При принятом здесь определении Early и Late
    // положительная ошибка означает, что локальный код
    // смещён в сторону поздней фазы.
    state.dll_integrator_chips_s -=dll_k2 * dll_phase_scale *state.dll_error * block_time;

    const double dll_rate_correction =-dll_k1 * dll_phase_scale * state.dll_error+ state.dll_integrator_chips_s;

    // 4. Параметры PLL.
    const double pll_wn = loop_wn(PLL_BANDWIDTH_HZ);
    const double pll_k1 =2.0 * LOOP_DAMPING * pll_wn;

    const double pll_k2 =pll_wn * pll_wn;
    // Интегратор хранит поправку к угловой частоте, рад/с.
    state.pll_integrator_rad_s +=pll_k2 * state.pll_error * block_time;
    const double pll_frequency_correction_hz =(pll_k1 * state.pll_error+ state.pll_integrator_rad_s) / TWO_PI;

    // 5. Обновляем скорости генераторов для следующего блока.
    state.code_rate_chips_s =state.code_rate_nominal_chips_s+ dll_rate_correction;
    state.carrier_doppler_hz =state.carrier_doppler_nominal_hz+ pll_frequency_correction_hz;
    return corr;
}