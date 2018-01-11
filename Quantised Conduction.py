# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re


#Constants of calibration curve
c = -3.7E-4
c_e = 0.3E-4
a = 5E-7
a_e = 3E-7
b = 3.7E-4
b_e = 0.2E-4
t1 = 0.49
t1_e = 0.03
t2 = 2.1
t2_e = 0.1

G0 = 7.7481E-5

def natural(text):
    # Function to define natural sorting
     return int(text) if text.isdigit() else text

def natural_sort(text):
    # Function to implement natural sorting
    return [ natural(c) for c in re.split('(\d+)', text) ]

def calibration(V):
    # Function defining the calibration equation to transform data from voltage
    # to quantised conductance
    G = a*np.exp(-V/t1)+b*np.exp(-V/t2)+c
    Out = G/G0
    return Out
    
def data_cal(data):
    # Function to implement the calibration on the data
    cal = np.genfromtxt(data, float, delimiter=',')
    cal = cal[(cal<=-0.39) & (cal>=-3.37)]
    cal_data = calibration(cal)
    return cal_data
  

def trim(filename, input_name, out_folder, note):
    # Function to trim the recorded data so that only traces of breaking events
    # based on specified inputs will be used. The resulting data are saved as 
    # text files corresponding to a single breaking event
    file = os.path.join(os.getcwd(), out_folder, filename)
    inp = os.path.join(os.getcwd(), out_folder, input_name)
    bound = np.genfromtxt(inp, dtype=int, delimiter = '\n')
    data = np.genfromtxt(file, delimiter = '\n')
    for i in range(0, len(bound)-1, 2):
        start = bound[i]
        end = bound[i+1]
        name = 'event' + str(note) + '.txt'
        output = os.path.join(os.getcwd(), out_folder, name)
        out = data[start:end]
        np.savetxt(output, out, newline = '\n') 
        note += 1
 
def produce_events(folder):
    # Function to isolate and save traces of breaking events from the recorded data
    file_list = glob.glob(os.path.join(os.getcwd(), folder, 'file*'))
    in_num = 1
    note = 1
    file_list.sort(key=natural_sort)
    for file in file_list:
        in_name = 'input' + str(in_num) + '.txt'
        trim(file, in_name, folder, note)
        file_len = glob.glob(os.path.join(os.getcwd(), folder, 'event*'))
        note = len(file_len) + 1
        in_num += 1
        
def convert(folder):
    # Function to convert the voltage traces of the produced breaking events
    # into conductance traces
    file_list = glob.glob(os.path.join(os.getcwd(), folder, 'event*'))
    count = len(file_list)
    print(count)
    c_files = []
    for file_path in file_list:
        c_files.append(data_cal(file_path))                                                               
    return c_files  
     
def combining_events(folder_name, name, lim):
    # Function to combine and save all of the breaking events into a single text file
    data = np.concatenate(convert(folder_name))
    data = data[(data<=lim)]
    data.sort()
    output = os.path.join(os.getcwd(), folder_name, name)
    np.savetxt(output, data, newline = '\n')
         

    
def plot(data_folder, bins, name):
    # Function to plot a histogram of the combined events
    file = glob.glob(os.path.join(os.getcwd(), data_folder, 'final*'))
    data = np.genfromtxt(file[0], delimiter = '\n')
    print(len(data))
    hist = plt.hist(data, bins, (1,10), color='#000000', histtype='barstacked')
    y_ax = hist[0]
    plt.xlabel("Conductance (G0)")
    plt.ylabel("Occurance")
    plt.xlim([1,10])
    plt.ylim([0,y_ax.max()])
    plt.xticks([1,2,3,4,5,6,7,8,9,10])
    xposition = [1,2,3,4,5,6,7,8,9,10]
    for xl in xposition:
        plt.axvline(x=xl, color='k', linestyle='--')
    output_name = os.path.join(os.getcwd(), data_folder, name)
    plt.savefig(output_name)
    
    
def plotevent(folder):
    # Function to plot the conductance traces of the individual breaking events
    file_list = glob.glob(os.path.join(os.getcwd(), folder, 'event*'))
    file_list.sort(key=natural_sort)
    count = 1
    for file in file_list:
        y = data_cal(file)
        y = y[(y<=5.0)]
        x = []
        if (y.size == 0):
            count += 1
        else:
            for i in range(len(y)):
                x.append(i/30)
            if not x:
                x.append(0)
            xarr = np.asarray(x)
            plt.figure(count)
            plt.xlim([0,np.amax(xarr)])
            plt.ylim([0,5])
            plt.plot(xarr,y)
            name = 'plotevent' + str(count) + '.png'
            output = os.path.join(os.getcwd(), folder, name)
            plt.savefig(output)
            plt.close()
            count += 1
        
def fourier(data_folder, bins_low, bins_up, step):
    # Function to perform fourier analysis on the binned data, calculate the signal's
    # power and its maximum peak. The height of the peak is then plotted against the 
    # bin size and the optimal size is returned. 
    file = glob.glob(os.path.join(os.getcwd(), data_folder, 'final*'))
    data = np.genfromtxt(file[0], delimiter = '\n')
    f_height = []
    bin_num = []
    n = bins_low
    while (n<bins_up):
        hist = plt.hist(data, n, (1,10), color='#000000', histtype='barstacked')
        plt.close()
        x_d = hist[1]
        y_d = hist[0]
        spacing = x_d[1]-x_d[0]
        f_data = np.fft.fft(y_d)
        freq= np.fft.fftfreq(len(y_d), spacing)
        p_data = np.abs(f_data)**2.0
        i = freq > 0
        g = 0
        index = []
        power = p_data[i]
        while (power[g+1] - power[g] < 0):
            index.append(g)
            g += 1
        power_fin = np.delete(power, index)
        f_height.append(10*np.log10(power_fin.max()))
        bin_num.append(n)
        n+=step
    plt.plot(bin_num, f_height)
    name = 'fourier ' + data_folder + '.txt'
    output = os.path.join(os.getcwd(), data_folder, name)
    np.savetxt(output, np.c_[bin_num, f_height], newline = '\n')
    ind = np.argmax(f_height)
    print(bin_num[ind])

    

    
        

    
    
