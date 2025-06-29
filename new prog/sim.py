import numpy as np
import matplotlib.pyplot as plt

# Параметры
symbol_duration = 20e-6  # 20 мкс
sampling_rate = 100e6   # 100 MHz
num_samples = int(symbol_duration * sampling_rate) # Отсчетов в такте
harmonics_range = range(1, 10)  # Увеличили число гармоник
num_symbols = 10000  # Увеличено количество символов еще больше для лучшей статистики при данном SNR
snr_db = 3      # SNR в dB, установлено на 3 dB
snr_linear = 10**(snr_db / 10)
carrier_frequency = sampling_rate / 4 # Несущая частота (25 MHz)

# 16-QAM созвездие (Нормализованное)
constellation = [(-3-3j)/np.sqrt(18), (-3-1j)/np.sqrt(18), (-3+1j)/np.sqrt(18), (-3+3j)/np.sqrt(18),
                 (-1-3j)/np.sqrt(18), (-1-1j)/np.sqrt(18), (-1+1j)/np.sqrt(18), (-1+3j)/np.sqrt(18),
                 ( 1-3j)/np.sqrt(18), ( 1-1j)/np.sqrt(18), ( 1+1j)/np.sqrt(18), ( 1+3j)/np.sqrt(18),
                 ( 3-3j)/np.sqrt(18), ( 3-1j)/np.sqrt(18), ( 3+1j)/np.sqrt(18), ( 3+3j)/np.sqrt(18)]

def generate_carrier(carrier_frequency, sampling_rate, num_samples):
    """Генерирует комплексную синусоидальную несущую."""
    t = np.arange(num_samples) / sampling_rate
    carrier = np.exp(2j * np.pi * carrier_frequency * t)
    return carrier

def generate_qam_signal(symbol, num_samples, carrier):
    """Генерирует 16-QAM сигнал для заданного символа путем модуляции несущей."""
    # Модулируем несущую комплексным значением символа
    signal = symbol * carrier
    return signal

def obpf(signal, num_harmonics, num_samples, sampling_rate, carrier_frequency):
    """Выполняет ОБПФ вокруг несущей с заданной полушириной в гармониках FFT."""
    fft_values = np.fft.fft(signal)
    freq = np.fft.fftfreq(num_samples, d=1/sampling_rate) # Получаем частоты для точек БПФ
    freq_step = sampling_rate / num_samples # Шаг частоты FFT

    # Создаем маску для полосового фильтра вокруг несущей и ее отрицательной частоты
    mask = np.zeros(num_samples, dtype=bool)

    # Определяем полосу пропускания вокруг положительной и отрицательной несущей
    # Полуширина полосы в Гц: num_harmonics * freq_step
    positive_band_center = carrier_frequency
    negative_band_center = -carrier_frequency
    bandwidth_half_hz = num_harmonics * freq_step

    # Включаем бины, находящиеся в полосе пропускания
    mask[(freq >= (positive_band_center - bandwidth_half_hz)) & (freq <= (positive_band_center + bandwidth_half_hz))] = True
    mask[(freq >= (negative_band_center - bandwidth_half_hz)) & (freq <= (negative_band_center + bandwidth_half_hz))] = True

    # Применяем маску к FFT
    filtered_fft = fft_values * mask

    # Выполняем обратное преобразование
    reconstructed_signal = np.fft.ifft(filtered_fft)

    return reconstructed_signal

def demodulate_signal(received_signal, constellation, num_samples, carrier):
    """Демодулирует сигнал, находя ближайшую точку в созвездии (Minimum Distance) используя корреляцию."""
    # Для корреляционной демодуляции, сравниваем принятый сигнал с ожидаемыми сигналами для каждого символа
    min_distance = float('inf')
    closest_symbol_index = -1

    # Генерируем ожидаемый сигнал для каждой точки созвездия (символ * несущая)
    expected_signals = [symbol * carrier for symbol in constellation]

    # Для каждого ожидаемого сигнала, вычисляем расстояние (или корреляцию) с принятым сигналом
    for i, expected_signal in enumerate(expected_signals):
        # Расстояние по MSE (Mean Squared Error) эквивалентно минимизации, связанной с корреляцией
        # Distance = ||received - expected||^2 = ||received||^2 - 2*Re(received * expected*) + ||expected||^2
        # Поскольку ||received||^2 одинаково для всех символов, минимизация сводится к максимизации Re(received * expected*)
        # Или, в данном случае, для MSE: Sum(|received[n] - expected[n]|^2)
        # Среднее значение MSE по всем отсчетам
        distance = np.mean(np.abs(received_signal - expected_signal)**2)

        if distance < min_distance:
            min_distance = distance
            closest_symbol_index = i

    return closest_symbol_index

def add_awgn(signal, snr_linear):
    """Добавляет аддитивный белый гауссовский шум (AWGN) к комплексному сигналу."""
    # Рассчитываем мощность сигнала как среднюю мощность по всем отсчетам
    signal_power = np.mean(np.abs(signal)**2)
    noise_power = signal_power / snr_linear
    # Генерируем комплексный гауссовский шум
    noise = np.sqrt(noise_power/2) * (np.random.randn(len(signal)) + 1j * np.random.randn(len(signal)))
    noisy_signal = signal + noise
    return noisy_signal

# Симуляция
ber_results = {}
carrier = generate_carrier(carrier_frequency, sampling_rate, num_samples) # Генерируем несущую один раз
for num_harmonics in harmonics_range:
    errors = 0
    for i in range(num_symbols):
        # 1. Случайно выбираем символ
        symbol_index = np.random.randint(0, 16)
        true_symbol = constellation[symbol_index]

        # 2. Генерируем сигнал (модулированный)
        signal = generate_qam_signal(true_symbol, num_samples, carrier) # Передаем несущую

        # 3. Добавляем шум
        noisy_signal = add_awgn(signal, snr_linear)

        # 4. Выполняем ОБПФ (полосовой вокруг несущей)
        reconstructed_signal = obpf(noisy_signal, num_harmonics, num_samples, sampling_rate, carrier_frequency) # Передаем параметры фильтра

        # 5. Демодулируем
        demodulated_index = demodulate_signal(reconstructed_signal, constellation, num_samples, carrier)

        # 6. Проверяем на ошибку
        if demodulated_index != symbol_index:
            errors += 1

    # 7. Вычисляем BER
    ber = errors / num_symbols
    ber_results[num_harmonics] = ber
    print(f"Harmonics: {num_harmonics}, BER: {ber:.8f}")