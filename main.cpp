#include <iostream>
#include <vector>
#include <complex>
#include <cstdint>
#include <array>
#include <cmath>

#include "math/iq_samples.h"
#include "math/math_func.h"
#include "acquistion/acquisition.h"
#include "tracking/tracking.h"

int main()
{
    // 1. Параметры сигнала
    constexpr double fs = 62314000.0;
    constexpr double intermediate_frequency = 31157000.0;

    constexpr std::size_t samples_per_ms =static_cast<std::size_t>(fs / 1000.0);

    // 2. Читаем входной сигнал
    IQSamples signal;
    signal.sampling_frequency_hz = fs;

    if (!signal.load("Glonass_L1_CA_SIM_Fs_62Msps_4ms.dat"))
    {
        std::cerr << "Cannot load IQ samples\n";
        return 1;
    }

    std::cout << "Loaded samples: "<< signal.size() << '\n';

    // 3. Переносим сигнал в основную полосу
    for (std::size_t n = 0; n < signal.samples.size(); ++n)
    {
        const double phase = -2.0 * M_PI * intermediate_frequency* static_cast<double>(n) / fs;
        const std::complex<float> carrier(static_cast<float>(std::cos(phase)),static_cast<float>(std::sin(phase)));

        signal.samples[n] *= carrier;
    }

    // 4. Генерируем локальный код для Acquisition
    auto local_code = generate_glonass_code(static_cast<int32_t>(fs), 0);

    // 5. Проверяем длину сигнала
    if (signal.size() < samples_per_ms)
    {
        std::cerr << "Not enough samples for acquisition\n";
        return 1;
    }

    // 6. Acquisition: первый миллисекундный блок
    std::vector<std::complex<float>> acquisition_block(signal.samples.begin(),signal.samples.begin() + samples_per_ms
    );

    AcquisitionResult acquisition_result = acquistion(
        acquisition_block,
        local_code,
        static_cast<float>(fs),
        -10000.0f,
        10000.0f,
        500.0f,
        0.005f
    );

    if (!acquisition_result.detected)
    {
        std::cerr << "Signal not detected. Tracking skipped.\n";
        return 1;
    }

    std::cout << "\nAcquisition successful\n";
    std::cout << "Doppler: "<< acquisition_result.doppler << " Hz\n";
    std::cout << "Delay: "<< acquisition_result.delay << " samples\n";

    // 7. Подготавливаем локальный код для Tracking
    // Tracking ожидает 511 чипов, а не отсчёты дискретизированного кода.
   const auto tracking_code_complex =generate_glonass_code(511000, 0);

std::array<int, 511> tracking_code{};

for (std::size_t i = 0; i < tracking_code.size(); ++i)
{
    tracking_code[i] =tracking_code_complex[i].real() >= 0.0f ? 1 : -1;
}
    // 8. Инициализируем состояние Tracking
    TrackingState tracking_state{};

    tracking_state.carrier_doppler_hz =acquisition_result.doppler;

    tracking_state.code_phase_chips = static_cast<double>(acquisition_result.delay)* 511000.0 / fs;

    // 9. Tracking: обрабатываем последующие миллисекундные блоки
    for (std::size_t start = samples_per_ms;start + samples_per_ms <= signal.size(); start += samples_per_ms)
    {
        std::vector<std::complex<float>> tracking_block(
            signal.samples.begin() + start,
            signal.samples.begin() + start + samples_per_ms
        );

        auto tracking_result = tracking_process_block(tracking_block,tracking_code,tracking_state,fs
        );

        std::cout << "\nTracking block at sample " << start << '\n';
        std::cout << "Prompt: " << tracking_result.prompt << '\n';
        std::cout << "DLL error: "<< tracking_state.dll_error << '\n';
        std::cout << "PLL error: "<< tracking_state.pll_error << '\n';
    }

    return 0;
}
