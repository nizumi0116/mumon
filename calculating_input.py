import numpy as np
from scipy.fftpack import fft, ifft, fftfreq
from scipy import interpolate
from matplotlib import pyplot as plt

ax = plt.gca()

data_test = 2050
data_over = 329

time_test = np.arange(1*10**(-9), data_test*4*10**(-9), 4*10**(-9))

data_array = np.load('ratio.npz')
ratio = data_array['x']
angle_dif = data_array['y']

dt = 4*10**(-9)
freq_fft_test = fftfreq(data_test, dt)
freq_fft_test_plus = freq_fft_test[0:data_test//2]

freq_array = np.arange(1*10**6, 41*10**6, 1*10**6)
freq_interpolate = np.arange(0, 41*10**6, 0.5*10**6)
func_amp = interpolate.InterpolatedUnivariateSpline(freq_array, ratio, k = 1)
func_angle = interpolate.InterpolatedUnivariateSpline(freq_array, angle_dif, k = 1)

result_amp_test = func_amp(freq_fft_test_plus)
for i in range(data_over, data_test//2):
    result_amp_test[i] = result_amp_test[data_over-1]

result_amp_test_min = result_amp_test[data_over-1]
result_angle_test = func_angle(freq_fft_test_plus)

for i in range(data_over, data_test//2):
    result_angle_test[i] = result_angle_test[data_over-1]

result_angle_test_min = result_angle_test[data_over-1]

run_array = [75, 78]

for run in run_array:
    if run == 83:
        path_test = f'rawdata/run_000{run}/run000{run}_ch1.txt'
    else:
        path_test = f'rawdata/run_000{run}/run000{run}_ch0.txt'
    
    l_test = np.loadtxt(path_test)
    if run == 73 or run == 76 or run == 79 or run == 83 or run == 82:
        l_test = l_test * 2.51
    if run == 50 or run == 51 or run == 52 or run == 54 or run == 59 or run == 66 or run == 69 or run == 72 or run == 74 or run == 75 or run == 77 or run == 78 or run == 81:
        l_test = l_test/10
    l_test = l_test.reshape(-1, data_test+2)
    l_test = l_test[:, 2:]
    
    #move baseline
    test_ave = np.mean(l_test[:, 0:200], axis = 1)
    for i in range(l_test.shape[0]):
        l_test[i:i+1, :] = l_test[i:i+1, :] - test_ave[i]
        
    #fft
    lf_test = fft(l_test)
    lf_test_real = lf_test.real
    nan_array = np.where(lf_test_real == 0)
    nan_array = np.array(nan_array)
    print(nan_array)
    lf_test_imag = lf_test.imag

    #delete nan array
    if nan_array.shape[1] != 0:
        row = 0
        for i in range(nan_array.shape[1]):
            print("delete: ", nan_array[0, i])
            lf_test = np.delete(lf_test, nan_array[0, i] - row, 0)
            lf_test_real = np.delete(lf_test_real, nan_array[0, i] - row, 0)
            lf_test_imag = np.delete(lf_test_imag, nan_array[0, i] - row, 0)
            row = row + 1
        print("size: ", lf_test.shape)

    print(np.where(lf_test_real == 0))
    
    #reproducing
    test_angle = np.arctan(lf_test_imag/lf_test_real)
    test_amp = lf_test_real/(np.cos(test_angle))

    test_amp_plus = test_amp[:, 0:data_test//2]
    test_angle_plus = test_angle[:, 0:data_test//2]

    for i in range(lf_test.shape[0]):
        if i % 100 is 0:
            print("processed:", i, "events")
        rep_test_complex = []
        rep_test_real = test_amp_plus[i:i+1, :] * 1/result_amp_test * np.cos(test_angle_plus[i:i+1, :] - result_angle_test)
        rep_test_imag = test_amp_plus[i:i+1, :] * 1/result_amp_test * np.sin(test_angle_plus[i:i+1, :] - result_angle_test)
        rep_test_real_min = test_amp[i:i+1, np.argmin(freq_fft_test):np.argmin(freq_fft_test)+1] * 1/result_amp_test_min * np.cos(test_angle[i:i+1, np.argmax(freq_fft_test):np.argmax(freq_fft_test)+1] - result_angle_test_min)
        rep_test_imag_min = -test_amp[i:i+1, np.argmin(freq_fft_test):np.argmin(freq_fft_test)+1] * 1/result_amp_test_min * np.sin(test_angle[i:i+1, np.argmax(freq_fft_test):np.argmax(freq_fft_test)+1] - result_angle_test_min)

        for x in range(0, data_test//2):
            rep_test_complex = np.append(rep_test_complex, complex(rep_test_real[:, x:x+1], rep_test_imag[:, x:x+1]))
        rep_test_complex = np.append(rep_test_complex, complex(rep_test_real_min, rep_test_imag_min))
        rep_test_real = np.append(rep_test_real, rep_test_real_min)
        rep_test_imag = np.append(rep_test_imag, rep_test_imag_min)

        for x in range(data_test//2-1, 0, -1):
            rep_test_real = np.append(rep_test_real, rep_test_real[x])
            rep_test_imag = np.append(rep_test_imag, -rep_test_imag[x])
            rep_test_complex= np.append(rep_test_complex, complex(rep_test_real[x], -rep_test_imag[x]))

        if i is 0:
            rep_test_comp = rep_test_complex
            rep_test_real_comp = rep_test_real
            rep_test_imag_comp = rep_test_imag

        else:
            rep_test_comp = np.vstack((rep_test_comp, rep_test_complex))
            rep_test_real_comp = np.vstack((rep_test_real_comp, rep_test_real))
            rep_test_imag_comp = np.vstack((rep_test_imag_comp, rep_test_imag))

    #try drawing
    rep_test_real_comp = rep_test_real_comp[150:151, :].reshape(data_test, -1)
    #plt.plot(freq_fft_test, rep_test_real_comp)
    #plt.show()

    rep_test = ifft(rep_test_comp)
    rep_test_real = rep_test.real
    np.savetxt(f'process/run_000{run}/run000{run}_ch0.dat', rep_test.real, delimiter = '\n', fmt = '%.18f')
    #print(rep_test.shape)
    l_test = l_test[150:151, :].reshape(data_test, -1)
    #plt.plot(time_test * 1000000000, l_test, color = "tab:green", label = 'output')
    rep_test = rep_test[150:151, :].reshape(data_test, -1)
    #plt.plot(time_test * 1000000000, rep_test, color = "tab:orange", label = 'calc input')
    #ax.legend(loc = 0)
    #plt.show()

    print("processed run", run)
