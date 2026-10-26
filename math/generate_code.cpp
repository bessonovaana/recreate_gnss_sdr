#include "math_func.h"
std::vector<std::complex<float>> generate_glonass_code(
    int32_t fs,
    uint32_t chip_shift)
{
    constexpr int L = 511;
    constexpr int32_t fc = 511000;


    // Число отсчётов на период кода — как в оригинале GNSS-SDR.
    const int32_t samples_per_code =static_cast<int32_t>(static_cast<double>(fs) /(static_cast<double>(fc) / L));

    auto R = std::bitset<9>{}.set();
    std::array<int, L> C{};

    // Генерация последовательности из 511 чипов.
    for (int n = 0; n < L; ++n){
        const bool outputbit = R[2];
        const bool feedback1 = R[4] ^ R[0];

        for (int k = 0; k < 8; ++k){
            R[k] = R[k + 1];
        }

        R[8] = feedback1;
        C[n] = outputbit ? 1 : -1;
    }

    // Циклический сдвиг последовательности.
    std::array<int, L> C_shift{};

    for (int n = 0; n < L; ++n){
        C_shift[n] = C[(n + chip_shift) % L];
    }

    std::vector<std::complex<float>> D(samples_per_code);

    // Дискретизация кода.
    for (int i = 0; i < samples_per_code; ++i)
    {
        const int chip_index = static_cast<int>(std::floor(static_cast<double>(i) * fc /static_cast<double>(fs)));
        // Защита от выхода за границы массива.
        const int index = std::min(chip_index, L - 1);

        D[i] = std::complex<float>(
            static_cast<float>(C_shift[index]), 0.0f);
    }

    // Как в оригинале: последний отсчёт относится к последнему чипу.
    D.back() = std::complex<float>(
        static_cast<float>(C_shift[L - 1]), 0.0f);

    return D;
}