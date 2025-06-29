#include <iostream>
#include <chrono>
#include <fftw3.h>
#include <cstring>
#include <math.h>
// #include <sys/resource.h>

int N;
double normalizeCoef;	// Нормализующий коэффициент
uint64_t COUNT_TEST;

void gener_start_mass(fftw_complex *symbolCoding)
{
    int i;
    for (i = 0; i < 16; ++i)
    switch (i)
    {
        case 0b0000: symbolCoding[i][0] = -3.0; symbolCoding[i][1] =  3.0; break; // 0
        case 0b0001: symbolCoding[i][0] = -1.0; symbolCoding[i][1] =  3.0; break; // 1
        case 0b0010: symbolCoding[i][0] = -3.0; symbolCoding[i][1] =  1.0; break; // 2
        case 0b0011: symbolCoding[i][0] = -1.0; symbolCoding[i][1] =  1.0; break; // 3
        case 0b0100: symbolCoding[i][0] = -3.0; symbolCoding[i][1] = -3.0; break; // 4
        case 0b0101: symbolCoding[i][0] = -1.0; symbolCoding[i][1] = -3.0; break; // 5
        case 0b0110: symbolCoding[i][0] = -3.0; symbolCoding[i][1] = -1.0; break; // 6
        case 0b0111: symbolCoding[i][0] = -1.0; symbolCoding[i][1] = -1.0; break; // 7
        case 0b1000: symbolCoding[i][0] =  3.0; symbolCoding[i][1] =  3.0; break; // 8
        case 0b1001: symbolCoding[i][0] =  1.0; symbolCoding[i][1] =  3.0; break; // 9
        case 0b1010: symbolCoding[i][0] =  3.0; symbolCoding[i][1] =  1.0; break; // 10
        case 0b1011: symbolCoding[i][0] =  1.0; symbolCoding[i][1] =  1.0; break; // 11
        case 0b1100: symbolCoding[i][0] =  3.0; symbolCoding[i][1] = -3.0; break; // 12
        case 0b1101: symbolCoding[i][0] =  1.0; symbolCoding[i][1] = -3.0; break; // 13
        case 0b1110: symbolCoding[i][0] =  3.0; symbolCoding[i][1] = -1.0; break; // 14
        case 0b1111: symbolCoding[i][0] =  1.0; symbolCoding[i][1] = -1.0; break; // 15
        default:   symbolCoding[i][0] = 0.0; symbolCoding[i][1] = 0.0; // Обработка неожиданной ситуации (опционально)
    }

    for (; i < N; ++i)
    {
        symbolCoding[i][0] = 0.0; symbolCoding[i][1] = 0.0;
    }
}



void OFFTW_func(fftw_plan planBack, fftw_complex* offtDataOut)
{
    int i;
    // ОБПФ
	fftw_execute(planBack);
	// Нормализация. В будущем заменить в зависимости от параметров приемника-передатчика
	for (i = 0; i < N; ++i)
	{
		offtDataOut[i][0] *= normalizeCoef;
		offtDataOut[i][1] *= normalizeCoef;
	}
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
    uint64_t k = 0;

    fftw_complex *offtDataIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *offtDataOut = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

    memset(offtDataIn, 0, sizeof(fftw_complex) * N);
    memset(offtDataOut, 0, sizeof(fftw_complex) * N);

    fftw_plan planBack;
    // Заранее создаем план для ОБПФ. 
    planBack = fftw_plan_dft_1d(N, offtDataIn, offtDataOut, FFTW_BACKWARD, FFTW_EXHAUSTIVE | FFTW_CONSERVE_MEMORY);

    gener_start_mass(offtDataIn);

    // }

    startFunc = std::chrono::system_clock::now();

    // функция обработки
    // {
    for (k = 0; k < COUNT_TEST; ++k)
        OFFTW_func(planBack, offtDataOut);

    // }    


    endFunc = std::chrono::system_clock::now();
    std::cout 
    // << "Время работы функции: " << std::chrono::duration<double>(endFunc - startFunc).count() 
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