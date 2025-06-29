#include <iostream>
#include <chrono>
#include <fftw3.h>
#include <cstring>
#include <math.h>
#include <vector>
#include <complex>
// #include <sys/resource.h>


int N;
double normalizeCoef;	// Нормализующий коэффициент
uint64_t COUNT_TEST;


std::complex<double> short_modulation(unsigned char cadrs)
{
 std::complex<double> symbolCoding;
 switch (cadrs) {
  case 0b0000: symbolCoding.real(-3.0); symbolCoding.imag(3.0); break; // 0
  case 0b0001: symbolCoding.real(-1.0); symbolCoding.imag(3.0); break; // 1
  case 0b0010: symbolCoding.real(-3.0); symbolCoding.imag(1.0); break; // 2
  case 0b0011: symbolCoding.real(-1.0); symbolCoding.imag(1.0); break; // 3
  case 0b0100: symbolCoding.real(-3.0); symbolCoding.imag(3.0); break; // 4
  case 0b0101: symbolCoding.real(-1.0); symbolCoding.imag(3.0); break; // 5
  case 0b0110: symbolCoding.real(-3.0); symbolCoding.imag(1.0); break; // 6
  case 0b0111: symbolCoding.real(-1.0); symbolCoding.imag(1.0); break; // 7
  case 0b1000: symbolCoding.real( 3.0); symbolCoding.imag(3.0); break; // 8
  case 0b1001: symbolCoding.real( 1.0); symbolCoding.imag(3.0); break; // 9
  case 0b1010: symbolCoding.real( 3.0); symbolCoding.imag(1.0); break; // 10
  case 0b1011: symbolCoding.real( 1.0); symbolCoding.imag(1.0); break; // 11
  case 0b1100: symbolCoding.real( 3.0); symbolCoding.imag(3.0); break; // 12
  case 0b1101: symbolCoding.real( 1.0); symbolCoding.imag(3.0); break; // 13
  case 0b1110: symbolCoding.real( 3.0); symbolCoding.imag(1.0); break; // 14
  case 0b1111: symbolCoding.real( 1.0); symbolCoding.imag(1.0); break; // 15
  default: symbolCoding.real(0.0); symbolCoding.imag(0.0);
 }
 return symbolCoding;
}



int main(int argc, char ** argv)
{
    std::chrono::system_clock::time_point startFunc, endFunc, startProgram;
    std::chrono::duration<double> srTime;
    // struct rusage dataStatus;
    startProgram = std::chrono::system_clock::now();
    N = atoi(argv[1]);
    normalizeCoef = 1/sqrt(N);
    COUNT_TEST = atoi(argv[2]);
    // Стартовые настройки.
    // {
    uint64_t k;
    int i;
    unsigned char cadrs[16];
    
    std::vector<std::vector<std::vector<std::complex<double>>>> offtDataIn(17, std::vector<std::vector<std::complex<double>>>(N, std::vector<std::complex<double>>(N)));
     std::vector<std::complex<double>> offtDataOut(N);
     // Инициализируем начальный массив данных:
     for (i = 0; i < 17; ++i)
     {
      for (unsigned int n = 0; n < N; ++n)
      {
       for (unsigned int k = 0; k < N; ++k)
       {
        offtDataIn[i][n][k] = short_modulation(i) * std::exp(std::complex<double>(0, 2 * M_PI * k * n / static_cast<double>(N)));
       }
      }
     }

     for (i = 0; i < 16; ++i)
        cadrs[i] = i;
    // }

    startFunc = std::chrono::system_clock::now();

    // функция обработки
    // {
    for (k = 0; k < COUNT_TEST; ++k)
    {

        for (i = 0; i < N; ++i)
            offtDataOut[i] = offtDataIn[cadrs[0]][i][0];
        for (unsigned int l = 0; l < N; ++l)
        {
            for (i = 1; i < 16; ++i)
            {
                offtDataOut[l] += offtDataIn[cadrs[i]][l][i];
            }
            for (; i < N; ++i)
            {
                offtDataOut[l] += offtDataIn[16][l][i];
            }
            offtDataOut[l] *= normalizeCoef;
        }
    }
    // }    


    endFunc = std::chrono::system_clock::now();
    std::cout
    //  << "Время работы функции: " << std::chrono::duration<double>(endFunc - startFunc).count() 
    // << "\nВремя работы программы: " << std::chrono::duration<double>(std::chrono::system_clock::now() - startProgram).count()
    << "Время одного преобразования в микросекундах: " << std::chrono::duration<double>(endFunc - startFunc).count() * 1000000 / COUNT_TEST << std::endl;

    // if (getrusage(RUSAGE_SELF, &dataStatus) == -1)
    // {
    //     std::cout << "Статусные данные не получены" << std::endl;
    // }
    // else
    // {
    //     std::cout << "Max RAM: " << dataStatus.ru_maxrss << " KiB\n"
    //     << "Разделяемая память: " << dataStatus.ru_ixrss << "\n"
    //     << "Неразделяемая память: " << dataStatus.ru_idrss << "\n"
    //     << "Неразделяемый стек: " << dataStatus.ru_isrss << std::endl;

    // }
    
        return 0;
}