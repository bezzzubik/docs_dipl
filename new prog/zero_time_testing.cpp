#include <iostream>
#include <chrono>
#include <fftw3.h>
#include <cstring>
#include <math.h>
#include <vector>
#include <complex>
#include <sys/resource.h>


int N = 64;
double normalizeCoef;	// Нормализующий коэффициент
int COUNT_TEST = 100000;


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
    // }

    startFunc = std::chrono::system_clock::now();

    // функция обработки
    // {
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