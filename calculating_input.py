import numpy as np
from scipy.fftpack import fft, ifft, fftfreq
from scipy import interpolate
from matplotlib import pyplot as plt

ax = plt.gca()

data = 1024
data_plus = 512
data_test = 2050
data_over = 329
data_out_over = 164
time = np.arange(1*10**(-9), data*4*10**(-9), 4*10**(-9))
time_test = np.arange(1*10**(-9), data_test*4*10**(-9), 4*10**(-9))

data_array = np.load('ratio.npz')
ratio = data_array['x']
angle_dif = data_array['y']

dt = 4*10**(-9)
freq_fft = fftfreq(data, dt)
freq_fft_test = fftfreq(data_test, dt)
freq_fft_test_plus = freq_fft_test[0:data_test//2]
freq_fft_int = freq_fft[0:data_plus]
freq_fft_minus = freq_fft[data_plus:]
freq_array = np.arange(1*10**6, 41*10**6, 1*10**6)
freq_interpolate = np.arange(0, 41*10**6, 0.5*10**6)
func_amp = interpolate.InterpolatedUnivariateSpline(freq_array, ratio, k = 1)
func_angle = interpolate.InterpolatedUnivariateSpline(freq_array, angle_dif, k = 1)
result_amp = func_amp(freq_fft_int)
for i in range(data_out_over, data_plus):
    result_amp[i] = result_amp[data_out_over-1]
result_amp_min = result_amp[data_out_over-1]
result_angle = func_angle(freq_fft_int)
for i in range(data_out_over, data_plus):
    result_angle[i] = result_angle[data_out_over-1]
result_angle_min = result_angle[data_out_over-1]

result_amp_test = func_amp(freq_fft_test_plus)
for i in range(data_over, data_test//2):
    result_amp_test[i] = result_amp_test[data_over-1]
result_amp_test_min = result_amp_test[data_over-1]
result_angle_test = func_angle(freq_fft_test_plus)
for i in range(data_over, data_test//2):
    result_angle_test[i] = result_angle_test[data_over-1]
result_angle_test_min = result_angle_test[data_over-1]

run_array = [51, 52, 54, 56]

path_in = f'rawdata/fourier/input_40mhz.txt'
path_out = f'rawdata/fourier/output_40mhz.txt'
path_in_test = f'rawdata/input_60ns_15ns.txt'
path_out_test = f'rawdata/output_60ns_15ns.txt'
for run in run_array:
    if run == 83:
        path_test = f'rawdata/run_000{run}/run000{run}_ch1.txt'
    else:
        path_test = f'rawdata/run_000{run}/run000{run}_ch0.txt'

    #move baseline
    l_in = np.loadtxt(path_in)
    in_top = np.max(l_in)
    in_bottom = np.min(l_in)
    in_ave = (in_top + in_bottom)/2
    l_in = l_in - in_ave
    
    l_test = np.loadtxt(path_test)
    if run == 73 or run == 76 or run == 79 or run == 83 or run == 82:
        l_test = l_test * 2.51
    if run == 51 or run == 52 or run == 54 or run == 81:
        l_test = l_test/10
    l_test = l_test.reshape(-1, data_test+2)
    l_test = l_test[:, 2:]
    
    #reshape
    l_in = l_in.reshape(-1, data)
    l_in = l_in[0]
    
    #move baseline
    l_out = np.loadtxt(path_out)
    out_top = np.max(l_out)
    out_bottom = np.min(l_out)
    out_ave = (out_top + out_bottom)/2
    l_out = l_out - out_ave
    
    test_ave = np.mean(l_test[:, 0:200], axis = 1)
    for i in range(l_test.shape[0]):
        l_test[i:i+1, :] = l_test[i:i+1, :] - test_ave[i]
        
    #reshape
    l_out = l_out.reshape(-1, data)
    l_out = l_out[0]
        
    #fft
    lf_in = fft(l_in)
    lf_in_real = lf_in.real
    lf_in_imag = lf_in.imag
    lf_out = fft(l_out)
    lf_out_real = lf_out.real
    lf_out_imag = lf_out.imag
    lf_test = fft(l_test)
    lf_test_real = lf_test.real
    #print(np.where(lf_test_real == 0))
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
    out_angle = np.arctan(lf_out_imag/lf_out_real)
    out_amp = lf_out_real/(np.cos(out_angle))
    test_angle = np.arctan(lf_test_imag/lf_test_real)
    test_amp = lf_test_real/(np.cos(test_angle))

    out_amp_plus = out_amp[0:data_plus]
    out_angle_plus = out_angle[0:data_plus]
    test_amp_plus = test_amp[:, 0:data_test//2]
    test_angle_plus = test_angle[:, 0:data_test//2]

    rep_real = out_amp_plus * 1/result_amp * np.cos(out_angle_plus - result_angle)
    rep_imag = out_amp_plus * 1/result_amp * np.sin(out_angle_plus - result_angle)

    rep_real_min = out_amp[data_plus] * 1/result_amp_min * np.cos(out_angle[data_plus] - result_angle_min)
    rep_imag_min = -out_amp[data_plus] * 1/result_amp_min * np.sin(out_angle[data_plus] - result_angle_min)
    rep_complex = []
    for i in range(0, data_plus):
        rep_complex.append(complex(rep_real[i], rep_imag[i]))
    rep_complex.append(complex(rep_real_min, rep_imag_min))
    rep_real = np.append(rep_real, rep_real_min)
    rep_imag = np.append(rep_imag, rep_imag_min)
    for i in range(data_plus-1, 0, -1):
        rep_real = np.append(rep_real, rep_real[i])
        rep_imag = np.append(rep_imag, -rep_imag[i])
        rep_complex.append(complex(rep_real[i], -rep_imag[i]))

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

    rep_test_real_comp = rep_test_real_comp[150:151, :].reshape(data_test, -1)
    #plt.plot(freq_fft_test, rep_test_real_comp)
    #plt.show()
    
    rep = ifft(rep_complex)
    #plt.plot(time, l_in, label = 'input')
    #plt.plot(time, l_out, label = 'output')
    #plt.plot(time, rep, label = 'calc input')
    #ax.legend(loc = 0)
    #plt.show()

    rep_test = ifft(rep_test_comp)
    rep_test_real = rep_test.real
    np.savetxt(f'process/run_000{run}/run000{run}_ch0.dat', rep_test.real, delimiter = '\n', fmt = '%.18f')
    #print(rep_test.shape)
    #l_test = l_test[150:151, :].reshape(data_test, -1)
    #plt.plot(time_test, l_test, label = 'output')
    #rep_test = rep_test[150:151, :].reshape(data_test, -1)
    #plt.plot(time_test, rep_test, label = 'calc input')
    #ax.legend(loc = 0)
    #plt.show()

    print("processed run", run)
    
#plt.plot(freq_fft_test, lf_test_real, label = 'output')
#plt.plot(freq_fft_test, rep_test_real, label = 'calc input')
#ax.legend(loc = 0)
#plt.show()

#plt.plot(freq_fft_test, lf_test_imag, label = 'output')
#plt.plot(freq_fft_test, rep_test_imag, label = 'calc input')
#ax.legend(loc = 0)
#plt.show()

#plt.plot(freq_fft, lf_in_real, label = 'input')
#plt.plot(freq_fft, lf_out_real, label = 'output')
#plt.plot(freq_fft, rep_real, label = 'calc input')
#ax.legend(loc = 0)
#plt.show()
#plt.plot(freq_fft, lf_in_imag, label = 'input')
#plt.plot(freq_fft, lf_out_imag, label = 'output')
#plt.plot(freq_fft, rep_imag, label = 'calc input')
#ax.legend(loc = 0)
#plt.show()
