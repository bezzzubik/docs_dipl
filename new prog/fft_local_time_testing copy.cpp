#include <iostream>
#include <chrono>
#include <cstring>
#include <vector>
#include <complex>
#include <cmath>
#include <sys/resource.h>


int N = 64;
double normalizeCoef;	// Нормализующий коэффициент
int COUNT_TEST = 100000;

constexpr double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::vector<Complex> ComplexVector;

// Прямое БПФ (для проверки)
void fft(ComplexVector& x) {
    const size_t N = x.size();
    if (N <= 1) return;

    ComplexVector even(N/2);
    ComplexVector odd(N/2);
    
    for (size_t i = 0; i < N/2; ++i) {
        even[i] = x[2*i];
        odd[i] = x[2*i + 1];
    }
    
    fft(even);
    fft(odd);
    
    for (size_t k = 0; k < N/2; ++k) {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k + N/2] = even[k] - t;
    }
}

// Обратное БПФ
void ifft(ComplexVector& x) {
    // Комплексное сопряжение входных данных
    for (auto& val : x) {
        val = std::conj(val);
    }
    
    // Выполняем прямое БПФ
    fft(x);
    
    // Комплексное сопряжение и нормализация результатов
    const size_t N = x.size();
    for (auto& val : x) {
        val = std::conj(val) / static_cast<double>(N);
    }
}


void gener_start_mass(ComplexVector& symbolCoding)
{
    int i;
    for (i = 0; i < N; ++i)
    switch (i)
    {
        case 0b0000: symbolCoding[i].real(-3.0); symbolCoding[i].imag(3.0); break; // 0
        case 0b0001: symbolCoding[i].real(-1.0); symbolCoding[i].imag(3.0); break; // 1
        case 0b0010: symbolCoding[i].real(-3.0); symbolCoding[i].imag(1.0); break; // 2
        case 0b0011: symbolCoding[i].real(-1.0); symbolCoding[i].imag(1.0); break; // 3
        case 0b0100: symbolCoding[i].real(-3.0); symbolCoding[i].imag(3.0); break; // 4
        case 0b0101: symbolCoding[i].real(-1.0); symbolCoding[i].imag(3.0); break; // 5
        case 0b0110: symbolCoding[i].real(-3.0); symbolCoding[i].imag(1.0); break; // 6
        case 0b0111: symbolCoding[i].real(-1.0); symbolCoding[i].imag(1.0); break; // 7
        case 0b1000: symbolCoding[i].real( 3.0); symbolCoding[i].imag(3.0); break; // 8
        case 0b1001: symbolCoding[i].real( 1.0); symbolCoding[i].imag(3.0); break; // 9
        case 0b1010: symbolCoding[i].real( 3.0); symbolCoding[i].imag(1.0); break; // 10
        case 0b1011: symbolCoding[i].real( 1.0); symbolCoding[i].imag(1.0); break; // 11
        case 0b1100: symbolCoding[i].real( 3.0); symbolCoding[i].imag(3.0); break; // 12
        case 0b1101: symbolCoding[i].real( 1.0); symbolCoding[i].imag(3.0); break; // 13
        case 0b1110: symbolCoding[i].real( 3.0); symbolCoding[i].imag(1.0); break; // 14
        case 0b1111: symbolCoding[i].real( 1.0); symbolCoding[i].imag(1.0); break; // 15
        default: symbolCoding[i].real(0.0); symbolCoding[i].imag(0.0);
    }
}



int main(int argc, char ** argv)
{
    std::chrono::system_clock::time_point startFunc, endFunc, startProgram;
    std::chrono::duration<double> srTime;
    struct rusage dataStatus;
    startProgram = std::chrono::system_clock::now();
    N = atoi(argv[1]);
    normalizeCoef = 1/sqrt(N);
    COUNT_TEST = atoi(argv[2]);
    // Стартовые настройки.
    // {
    int k = 0;

    ComplexVector signal(N);
    gener_start_mass(signal);

    // }

    startFunc = std::chrono::system_clock::now();

    // функция обработки
    // {
    for (k = 0; k < COUNT_TEST; ++k)
        ifft(signal);

    // }    

    endFunc = std::chrono::system_clock::now();
    std::cout << "Время работы функции: " << std::chrono::duration<double>(endFunc - startFunc).count() 
    << "\nВремя работы программы: " << std::chrono::duration<double>(std::chrono::system_clock::now() - startProgram).count()
    << "\nВремя одного преобразования в микросекундах: " << std::chrono::duration<double>(endFunc - startFunc).count() * 1000000 / COUNT_TEST << std::endl;

        if (getrusage(RUSAGE_SELF, &dataStatus) == -1)
    {
        std::cout << "Статусные данные не получены" << std::endl;
    }
    else
    {
        std::cout << "Max RAM: " << dataStatus.ru_maxrss << " KiB\n"
        << "Разделяемая память: " << dataStatus.ru_ixrss << "\n"
        << "Неразделяемая память: " << dataStatus.ru_idrss << "\n"
        << "Неразделяемый стек: " << dataStatus.ru_isrss << std::endl;

    }
    return 0;
}