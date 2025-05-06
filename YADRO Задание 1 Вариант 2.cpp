#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cassert>
#include <random>
using namespace std;

class FFT {
public:
    // Прямое преобразование Фурье (БПФ)
    static void forward(vector<complex<double>>& x) {
        int n = x.size();
        if (n <= 1) return;
        // Делим массив на четные и нечетные элементы
        vector<complex<double>> even(n / 2);
        vector<complex<double>> odd(n / 2);
        for (int i = 0; i < n / 2; i++) {
            even[i] = x[i * 2];
            odd[i] = x[i * 2 + 1];
        }
        // Рекурсивные вызовы для даже и нечетных частей
        forward(even);
        forward(odd);
        // Комбинируем результаты
        for (int k = 0; k < n / 2; k++) {
            complex<double> t = polar(1.0, -2 * M_PI * k / n) * odd[k];
            x[k] = even[k] + t;
            x[k + n / 2] = even[k] - t;
        }
    }
    // Обратное преобразование Фурье
    static void inverse(vector<complex<double>>& x) {
        int n = x.size();
        for (auto& c : x) c = conj(c); // Конъюгируем
        forward(x); // Прямое БПФ
        for (auto& c : x) c = conj(c) / static_cast<double>(n); // Обратный масштаб
    }
};
int main() {
    const int N = 12;
    vector<complex<double>> input(N);

    mt19937 generator(42);
    uniform_real_distribution<double> distribution(-10.0, 10.0);

    for (int i = 0; i < N; ++i) {
        input[i] = complex<double>(distribution(generator), distribution(generator));
    }
    cout << "Input data:" << endl;
    for (const auto& c : input) {
        cout << c << endl;
    }

    vector<complex<double>> fft_output = input;
    FFT::forward(fft_output);
    cout << "\nFFT output:" << endl;
    for (const auto& c : fft_output) {
        cout << c << endl;
    }

    vector<complex<double>> ifft_output = fft_output;
    FFT::inverse(ifft_output);
    cout << "\nInverse FFT output:" << endl;
    for (const auto& c : ifft_output) {
        cout << c << endl;
    }

    cout << "\nError comparison:" << endl;
    for (size_t i = 0; i < input.size(); ++i) {
        double error = abs(input[i] - ifft_output[i]);
        cout << "Error for index " << i << ": " << error << endl;
    }
    return 0;
}
