COUNT_TEST=50000000
g++ fftw_time_testing.cpp -std=c++20 -lm -lfftw3 -O0 -o fftw.out
g++ fft_local_time_testing.cpp -std=c++20 -lm -O0 -o fft.out
g++ dpf_time_testing.cpp -std=c++20 -lm -O0 -o dpf.out

for param in 8 16 20 32 40 64; do
    echo "Count garmon = $param, Count tests = $COUNT_TEST"
    # echo "ZERO RESOURSE:"
    # ./a.out $param $COUNT_TEST
    # echo ""

    echo "Optimise DFT..."
    ./dpf.out $param $COUNT_TEST

    echo ""
    echo "standart FFT..."
    ./fft.out $param $COUNT_TEST

    echo ""
    echo "Lib FFT..."
    ./fftw.out $param $COUNT_TEST
done