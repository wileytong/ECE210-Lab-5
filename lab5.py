from pylab import *
from scipy.signal import butter, lfilter, freqz,firwin
from scipy.io import wavfile


fs = 44100.0 # Sampling rate used in computer, Hz
fs_high = 44100*10 # high sampling rate 
T = 0.1

nyq_rate = fs/2

time_display = 2 # ms 

def create_input(f):
    test_signal_high = create_sinusoid(fs_high,T,f) 
    return test_signal_high

def sampling(test_signal):
    signal_sampled = test_signal[::10]
    return signal_sampled

def create_sinusoid(fs,T,f):
    t = linspace(0,T,T*fs) 
    #t = linspace(0,0.1,0.1*fs) 
    test_signal = cos(2*pi*f*t)*0.2
    return test_signal

def system_output(test_signal,Processing ="No processing",impulse_response=0):
    if Processing == "No processing":
        output = test_signal
    else :
        if Processing == "IF Filter":
            
            #L = len(test_signal)
            #fft_taps = fft(impulse_response,L)
            #fft_y = fft(test_signal,L)
            #fft_output =  fft_taps*fft_y
            
            #output = ifft(fft_output,L)
            
            output = lfilter(impulse_response,1,test_signal)
        else : 
            if Processing == "IF Filter+Envelope Detector":
                filtered_signal = lfilter(impulse_response,1,test_signal)
                y_envelope = Envelope_Detector(filtered_signal)
                output =  y_envelope   
                
            else:
                if Processing == "IF Filter+Coherent Detector":
                    filtered_signal = lfilter(impulse_response,1,test_signal)
                    y_coherent = Coherent_Detector(filtered_signal)
                    #y_coherent = Coherent_Detector(test_signal)
                    output =  y_coherent
                    
                else:

                    print("No command found")
                    return

    return output


def create_bandpass_filter(fc_low,fc_high):
    cutoff_hz = array([fc_low,fc_high])
    n = 127
    taps = firwin(n,cutoff_hz/nyq_rate,pass_zero=False)
    return taps

def low_pass_filter(cut_off):
    # return impulse response of a low pass filter with given cut off frequency
    fc_low = cut_off
    cutoff_hz = array([fc_low])
    n = 127
    taps_lpf = firwin(n,cutoff_hz/nyq_rate)
    return taps_lpf

def Envelope_Detector(f):
    # return output of an envelope dector with given input
    f_abs = abs(f)
    cut_off = 4e3
    taps_lpf = low_pass_filter(cut_off)
    y_envelope = convolve(f_abs,taps_lpf)
    return y_envelope

def Coherent_Detector(f):
    # return output of an envelope dector with given input
    L = len(f)
    t = linspace(0,L/fs,L)
    fft_len = 2**16
    
    if L > fft_len:
        f_fft = abs(fft(f,fft_len))
        carrier_freq = argmax(f_fft[:int(fft_len/2)])/fft_len*fs
    else :
        
        f_fft = abs(fft(f,L))
        carrier_freq = argmax(f_fft[:int(L/2)])/L*fs
    #carrier_freq = 14e3
    print("carrier_freq=%f Hz",carrier_freq)
    
    f_demod = f*sin(carrier_freq*2*pi*t)
    cut_off = 4e3
    taps_lpf = low_pass_filter(cut_off)
    y_coherent = convolve(f_demod,taps_lpf)
    return y_coherent


# Display functions 
def input_display(test_signal):
    npts_display = round(time_display*fs_high/1000)    
    t_axis = linspace(0,time_display,npts_display)
    plot(t_axis,test_signal[:npts_display])
    xlabel("Time / ms")
    ylabel("Amplitude")
    title("Input Signal $f(t)$ vs Time (ms)")
    
def time_domain_display_input(test_signal):
    test_signal = roll(test_signal,200)
    npts_display = round(time_display*fs/1000)    
    t_axis = linspace(0,time_display,npts_display)
    stem(t_axis,test_signal[:npts_display])
    xlabel("Time / ms")
    ylabel("Amplitude")
    title("Sampled Signal $f(nT)$ vs Time (ms)")

def time_domain_display_output(test_signal,mode="normal"):
    test_signal = roll(test_signal,200)
    L = len(test_signal)

    if mode=="fine tuned":
        y_fft = fft(test_signal)
        
        
        mono_freq = argmax(abs(y_fft[:int(L/2)]))/L*fs
        output = create_sinusoid(fs_high,T,mono_freq)
       
        npts_display = round(time_display*fs_high/1000)    
        t_axis = linspace(0,time_display,npts_display)
        plot(t_axis,output[:npts_display])
        
    if mode=="normal":
        npts_display = round(time_display*fs/1000)    
        t_axis = linspace(0,time_display,npts_display)
        plot(t_axis,test_signal[:npts_display])    

    xlabel("Time / ms")
    ylabel("Amplitude")
    title("Output Signal $y(t)$ vs Time (ms)")
    
def freq_domain_display(test_signal,unit="rms"):
    
    npts_display = 256
    
    y_fft = fft(test_signal,npts_display*2)
    
    freq_display = fs/2/1000
    freq_axis = linspace(0,freq_display,npts_display)
    if unit == "rms":
        plot(freq_axis,abs(y_fft)[:npts_display])
        ylabel("Magnitude")
    else :
        plot(freq_axis,20*log10(abs(y_fft)[:npts_display]))
        ylabel("Magnitude / dB")
    xlabel("Frequency / kHz")
    title("Magnitude vs Frequency (kHz)")
    #ylim([-30,50])
    
def freq_response(impulse_response):

    w, h = freqz(impulse_response, worN=8000)
    plot((w/pi)*nyq_rate/1000, 20*log10(abs(h)), linewidth=2)
    xlabel("Frequency / kHz")
    ylabel("Magnitude / dB")
    title("Filter frequency response")    
    ylim([-150,50])
    
    
    