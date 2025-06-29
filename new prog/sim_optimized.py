import numpy as np
import matplotlib.pyplot as plt

# Параметры
symbol_duration = 20e-6  # 20 мкс
sampling_rate = 100e6   # 100 MHz
num_samples = 200 # Отсчетов в такте (уменьшено до 200 для ускорения и изменения отношения полос)
harmonics_range = range(1, 129)  # Исследуем диапазон гармоник до 128
num_symbols = 10000  # Уменьшено число символов для ускорения
snr_db = -5      # SNR в dB, установлено на -5 dB для гарантии ошибок
snr_linear = 10**(snr_db / 10)
carrier_frequency = sampling_rate / 4 # Несущая частота (25 MHz)

# 16-QAM созвездие (Нормализованное)
constellation = [(-3-3j)/np.sqrt(18), (-3-1j)/np.sqrt(18), (-3+1j)/np.sqrt(18), (-3+3j)/np.sqrt(18),
                 (-1-3j)/np.sqrt(18), (-1-1j)/np.sqrt(18), (-1+1j)/np.sqrt(18), (-1+3j)/np.sqrt(18),
                 ( 1-3j)/np.sqrt(18), ( 1-1j)/np.sqrt(18), ( 1+1j)/np.sqrt(18), ( 1+3j)/np.sqrt(18),
                 ( 3-3j)/np.sqrt(18), ( 3-1j)/np.sqrt(18), ( 3+1j)/np.sqrt(16), ( 3+3j)/np.sqrt(18)] # Исправлена опечатка в последнем символе

def generate_carrier(carrier_frequency, sampling_rate, num_samples):
    """Генерирует комплексную синусоидальную несущую."""
    t = np.arange(num_samples) / sampling_rate
    carrier = np.exp(2j * np.pi * carrier_frequency * t)
    return carrier

def generate_qam_signal(symbol, carrier):
    """Генерирует 16-QAM сигнал для заданного символа путем модуляции несущей."""
    # Модулируем несущую комплексным значением символа
    signal = symbol * carrier
    return signal

def obpf(signal, num_harmonics, sampling_rate, carrier_frequency):
    """Выполняет ОБПФ вокруг несущей с заданной полушириной в гармониках FFT."""
    num_samples = len(signal)
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
    # Учитываем "перенос" отрицательных частот в верхнюю половину массива FFT
    positive_mask = (freq >= (positive_band_center - bandwidth_half_hz)) & (freq <= (positive_band_center + bandwidth_half_hz))
    negative_mask = (freq >= (negative_band_center - bandwidth_half_hz)) & (freq <= (negative_band_center + bandwidth_half_hz))

    mask = positive_mask | negative_mask

    # Применяем маску к FFT
    filtered_fft = fft_values * mask

    # Выполняем обратное преобразование
    reconstructed_signal = np.fft.ifft(filtered_fft)

    return reconstructed_signal

def demodulate_signal(received_signal, constellation, carrier):
    """Демодулирует сигнал, находя ближайшую точку в созвездии (Minimum Distance) используя корреляцию."""
    # Для корреляционной демодуляции, сравниваем принятый сигнал с ожидаемыми сигналами для каждого символа
    min_distance = float('inf')
    closest_symbol_index = -1

    # Генерируем ожидаемый сигнал для каждой точки созвездия (символ * несущая)
    expected_signals = [symbol * carrier for symbol in constellation]

    # Для каждого ожидаемого сигнала, вычисляем расстояние (или корреляцию) с принятым сигналом
    for i, expected_signal in enumerate(expected_signals):
        # Расстояние по MSE (Mean Squared Error) эквивалентно минимизации, связанной с корреляцией
        # Среднее значение MSE по всем отсчетам
        distance = np.mean(np.abs(received_signal - expected_signal)**2)

        if distance < min_distance:
            min_distance = distance
            closest_symbol_index = i

    return closest_symbol_index

def add_awgn(signal, snr_linear):
    """Добавляет аддитивный белый гауссовский шум (AWGN) к комплексному сигналу."""
    # Рассчитываем мощность сигнала как среднюю мощность по всем отсчетам
    # Поскольку созвездие нормализовано к средней мощности 1, и несущая имеет амплитуду 1,
    # средняя мощность символа в созвездии * несущая = 1 * 1^2 = 1.
    # Используем 1 как среднюю мощность символа для расчета мощности шума.
    signal_power = 1.0 # Средняя мощность символа после нормализации созвездия
    noise_power = signal_power / snr_linear
    # Генерируем комплексный гауссовский шум
    noise = np.sqrt(noise_power/2) * (np.random.randn(len(signal)) + 1j * np.random.randn(len(signal)))
    noisy_signal = signal + noise
    return noisy_signal

# Симуляция
ber_results = {}
carrier = generate_carrier(carrier_frequency, sampling_rate, num_samples) # Генерируем несущую один раз на весь такт символа

print("Начало симуляции...")

for num_harmonics in harmonics_range:
    errors = 0
    # Используем ту же несущую для всех символов внутри цикла гармоник
    for i in range(num_symbols):
        # 1. Случайно выбираем символ
        symbol_index = np.random.randint(0, 16)
        true_symbol = constellation[symbol_index]

        # 2. Генерируем сигнал (модулированный)
        signal = generate_qam_signal(true_symbol, carrier)

        # 3. Добавляем шум
        noisy_signal = add_awgn(signal, snr_linear)

        # 4. Выполняем ОБПФ (полосовой вокруг несущей)
        reconstructed_signal = obpf(noisy_signal, num_harmonics, sampling_rate, carrier_frequency)

        # 5. Демодулируем
        demodulated_index = demodulate_signal(reconstructed_signal, constellation, carrier)

        # 6. Проверяем на ошибку
        if demodulated_index != symbol_index:
            errors += 1

    # 7. Вычисляем BER
    ber = errors / num_symbols
    ber_results[num_harmonics] = ber
    print(f"Harmonics: {num_harmonics}, BER: {ber:.8f}")

print("Симуляция завершена.")

# Определение оптимального числа гармоник (минимальный BER)
min_ber = float('inf')
optimal_harmonics = -1

for harmonics, ber in ber_results.items():
    if ber < min_ber:
        min_ber = ber
        optimal_harmonics = harmonics

# Если минимальный BER > 0, то это найденный минимум.
# Если минимальный BER == 0, то либо 0 ошибок при всех гармониках (маловероятно при 1dB), либо
# минимум достигается при самой узкой полосе, которая дала 0 ошибок.
# Если все BER равны 0, возможно, SNR все еще слишком высок или число символов мало.

print(f"\nОптимальное число гармоник для минимального BER: {optimal_harmonics} (BER: {min_ber:.8f})")

# Визуализация результатов (опционально)
plt.figure()
plt.plot(list(ber_results.keys()), list(ber_results.values()), marker='o')
plt.xlabel('Число гармоник ОБПФ')
plt.ylabel('BER')
plt.title('Зависимость BER от числа гармоник ОБПФ')
plt.grid(True)
plt.show() 