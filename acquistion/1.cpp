// #include <gnuradio/io_signature.h>

#include <array>
#include <bitset>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <vector>
#include <complex>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "iq_samples.h"
#include "math_func.h"

struct AcquisitionResult
{
    float doppler;
    long long delay;
    float glrt;
    bool detected;
};







float calculate_power(const std::vector<std::complex<float>>& iq_samples)
{
    if (iq_samples.empty()){
        throw std::invalid_argument("Input signal is empty");
    }

    double sum = 0.0;
    for (std::size_t k = 0; k < iq_samples.size(); ++k){
        sum += std::norm(iq_samples[k]);
    }
    return static_cast<float>(sum / static_cast<double>(iq_samples.size()));
}




std::vector<std::complex<float>> doppler_shift (float f_d, const std::vector<std::complex<float>>& iqsamples, float fs){
    std::vector<std::complex<float>> x_d(iqsamples.size());
    float T_s = 1/fs;
  
    for (int k=0; k<iqsamples.size(); k++){
            std::complex<float> carrier = std::exp(std::complex<float>(0.0f, -2.0f*M_PI*f_d*k *T_s));
            x_d[k]=iqsamples[k]*carrier;
    }
    // std::cout<<"f_d = "<<f_d<<std::endl;
    return x_d;

}

std::vector<std::complex<float>> fft (const std::vector<std::complex<float>> &iq_samples_time){
    int N = 2;
    while (N<iq_samples_time.size()){
        N<<=1;
    }
    std::vector<std::complex<float>> sample_freq(N, {0.0f, 0.0f});
    std::copy(iq_samples_time.begin(),iq_samples_time.end(),sample_freq.begin());
    for (int i=1, j=0; i<N; i++){
        std::size_t bit =N>>1;
        while (j&bit){
            j^=bit;
            bit>>=1;
        }
        j^=bit;
        if (i<j){
            std::swap(sample_freq[i], sample_freq[j]);
        }
    }
    for (std::size_t len=2; len<=N; len<<=1){
        float angle = -2.0f*M_PI/static_cast<float>(len);
        std::complex<float> Wlen = std::polar(1.0f,angle);
        for (std::size_t i=0; i<N; i+=len){
            std::complex<float> W(1.0f,0.0f);
            for (std::size_t j=0;j<len/2;j++){
                std::complex<float> u=sample_freq[i+j];
                std::complex<float> v=sample_freq[i+j+len/2]*W;
                sample_freq[i+j]=u+v;
                sample_freq[i+j+len/2]=u-v;
                W*=Wlen;
            }
        }
    }
    return sample_freq;
}

std::vector<std::complex<float>> conjugate(const std::vector<std::complex<float>>& spectrum)
{
    std::vector<std::complex<float>> result(spectrum.size());
    for (std::size_t k = 0; k < spectrum.size(); k++)
    {
        result[k] = std::conj(spectrum[k]);
    }
    return result;
}

std::vector<std::complex<float>> multiply_spectra(const std::vector<std::complex<float>>& signal_fft,const std::vector<std::complex<float>>& code_fft)
{
    const std::size_t N = signal_fft.size();
    std::vector<std::complex<float>> result(N);
    for (std::size_t k = 0; k < N; k++){
        result[k] =signal_fft[k] *std::conj(code_fft[k]);
    }
    return result;
}

std::vector<std::complex<float>> ifft(const std::vector<std::complex<float>>& spectrum){
    const std::size_t N = spectrum.size();
    std::vector<std::complex<float>> conjugated(N);
    for (std::size_t k = 0; k < N; k++){
        conjugated[k] = std::conj(spectrum[k]);
    }
    std::vector<std::complex<float>> result =fft(conjugated);
    for (std::size_t n = 0; n < N; n++){
        result[n] =std::conj(result[n]) /static_cast<float>(N);
    }
    return result;
}



AcquisitionResult acquistion (
    const std::vector<std::complex<float>>& iq_samples,
    const std::vector<std::complex<float>>& local_code,
    float fs,
    float f_min,
    float f_max,
    float f_step,
    float threshold){

    float global_peak = 0.0f;
    float best_doppler = 0.0f;
    long long best_delay = 0;

    auto code_fft = fft(local_code);
    std::cout << "Input size = " << iq_samples.size() << '\n';
    std::cout << "Local code size = " << local_code.size() << '\n';
    std::cout << "Code FFT size = " << code_fft.size() << '\n';

    auto power_in = calculate_power(iq_samples);

    std::vector<std::vector<float>> power_grid;
    std::vector<float> doppler_values;
    for (float f_d = f_min; f_d <= f_max; f_d += f_step){
    auto x_d = doppler_shift(f_d, iq_samples, fs);
    auto signal_fft = fft(x_d);
    auto corr_spectrum = multiply_spectra(signal_fft, code_fft);
    auto corr = ifft(corr_spectrum);

    std::vector<float> power = correlation_power(corr);

    power_grid.push_back(power);
    doppler_values.push_back(f_d);

    // Пока ищем максимум только среди первых iq_samples.size()
    // отсчётов корреляции.
    const std::size_t valid_size =
        std::min(iq_samples.size(), power.size());

    auto max_it = std::max_element(power.begin(),power.begin() + valid_size);
    const std::size_t delay =std::distance(power.begin(), max_it);
    const float peak = *max_it;

    if (f_d == -9500.0f){
        std::cout << "CHECK: f_d = " << f_d<< ", peak = " << peak<< ", delay = " << delay << '\n';
    }

    if (peak > global_peak){
        global_peak = peak;
        best_doppler = f_d;
        best_delay = static_cast<long long>(delay);
    }
}

  
    std::cout<<"f_d = "<<best_doppler<<std::endl;
    std::cout<<"best delay = "<<best_delay<<std::endl;
    std::cout<<"best peak = "<<global_peak<<std::endl;


    const float glrt =2.0f * global_peak /(static_cast<float>(iq_samples.size()) * power_in);
    std::cout<<"gltr = "<<glrt<<std::endl;
    
    if (glrt > threshold)
    {
        std::cout << "SIGNAL DETECTED\n";
    }
    else
    {
        std::cout << "SIGNAL NOT DETECTED\n";
    }
    return {
        best_doppler,
        best_delay,
        glrt,
        glrt > threshold,
    };
}




// int main(){
//     IQSamples signal;
//     signal.sampling_frequency_hz = fs;
//     signal.load("Glonass_L1_CA_SIM_Fs_62Msps_4ms.dat");
//     auto local_code = generate_glonass_code(fs, 0);


//     const std::size_t samples_per_ms = fs / 1000;

//     std::vector<std::complex<float>> test_block(signal.samples.begin(),signal.samples.begin() + samples_per_ms);    
//     const float intermediate_frequency = 9540000.0f;

//     for (std::size_t n = 0; n < test_block.size(); ++n)
//     {
        
//         const float phase =-2.0f * static_cast<float>(M_PI) *intermediate_frequency * static_cast<float>(n) / fs;
//         test_block[n] *= std::complex<float>(std::cos(phase),std::sin(phase));
        
//     }

// std::cout << "Test block size = "<< test_block.size()<< std::endl;
//     AcquisitionResult res = acquistion(test_block,local_code,62314000.0f,-5000.0f,5000.0f,500.0f,0.005f);
//     return 0;
// }