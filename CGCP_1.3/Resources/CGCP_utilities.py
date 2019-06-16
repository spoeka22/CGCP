import time
import os

import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erf
from scipy.optimize import curve_fit
import pandas as pd

"""
Utility Tools
"""
#easy extraction tool. generates a vector with retentiontime and counts
def extract_data(sequence,injection_number, detector):
    #Convinient function to get at the data
    chfile = sequence.injections[injection_number].raw_files[detector]
    return chfile.times, chfile.values 

#extracts information given in the injection filename and info
def extract_sample_info(sequence):
    # this  should also find pressure from gas data to automatically include
    # that in the calibration
    # global concentration_list, name_list, calibration, injection, conc
    global injection
    concentration_list = []
    name_list = []
    calibration = False
    pressure_list = []
    for injection in sequence.injections:
        sample_info = injection.metadata['sample_info'][1:-5]
        sample_name_chromatogram = injection.metadata['sample_name']
        print("Sample info= " + str(sample_info))
        print("Sample name in chromatogram = " + str(sample_name_chromatogram))
        #find the right sample name to use for saving the plots/data

        if "µM" in sample_info:
            print("This is a calibration sequence.")
            calibration = True
            conc_info = sample_info[:sample_info.find("µM")]
            conc = int(''.join(filter(str.isdigit,
                                      conc_info)))  # this only works if there is only ONE number in the string before "µM" and that is the concentration
            concentration_list.append(conc)
        elif "AW_" in sample_info:
            print("AW in sample info - taken as sample name.")
            sample_name = sample_info[sample_info.find("AW_"):]
            name_list.append(sample_name)
        elif "Pd" in sample_name_chromatogram or "Au" in sample_name_chromatogram:  # the idea is that even if the naming if AW is not given, that the sample name is found
            print("Au or Pd sample. Name taken from chromatogram sample name.")
            sample_name = sample_name_chromatogram
            name_list.append(sample_name)
        elif "blank" in sample_name_chromatogram:
            print("blank sample  - name taken from chromatogram sample name.")
            sample_name = sample_name_chromatogram
            name_list.append(sample_name)
        else:  # not sure if that actually makes sense, but leave it for now (potentially needs to be changed, also depends on how systematically I name the samples
            print("Sample info chosen as name.")
            sample_name = sample_info
            name_list.append(sample_name)
        print("Selected sample name for evaluation: " + str(sample_name))
        #additional information: pressure (for gas injections)
        if "p=" in sample_name_chromatogram or "p =" in sample_name_chromatogram:
            pressure_string = sample_name_chromatogram
        elif "p=" in sample_info or "p =" in sample_info:
            pressure_string = sample_info
        else:
            pressure_string = "=0"
            # print(pressure_string)
            # print(pressure_string[pressure_string.find("=")+1:])
        if "mbar" in pressure_string:
            pressure = float(pressure_string[pressure_string.find("=")+1:pressure_string.find("mbar")])
        else:
            pressure = float(pressure_string[pressure_string.find("=")+1:])
        print(pressure)
        injection.metadata['pressure'] = pressure
        pressure_list.append(pressure)

    return concentration_list, name_list, calibration, pressure_list

#general function to save plots as png and/or pdf
def save_a_plot(output_path, plotname, pdf = True, png = True, dots_per_inch=300):
    try:
        if pdf == True:
            plt.savefig(output_path + plotname + '.pdf', dpi=dots_per_inch, bbox_inches='tight')
        if png == True:
            plt.savefig(output_path + plotname + '.png', dpi= dots_per_inch, bbox_inches='tight')
        plt.close()
    except OSError as err:
        print("OS error: {0}".format(err))
        print("Output path contains invalid characters.")
        new_filename = input("Enter new filename:")
        if pdf == True:
            plt.savefig(output_path + new_filename + '.pdf', dpi=dots_per_inch, bbox_inches='tight')
        if png == True:
            plt.savefig(output_path + new_filename + '.png', dpi= dots_per_inch, bbox_inches='tight')
        plt.close()

#step function
def errorfunc(time, center, size, height):
    time = np.array(time)
    step_speed = 10000 #making the erf-function a almost step-function
    return height+size*erf(step_speed*(time-center))

#curvefit erf function
def fit_erf_background(time_data, func_data, center):
    time_data=np.array(time_data)
    func_data=np.array(func_data)
    def erf_background(time, size, height): #the errorfunction is re-defined for parsing it correctly into "curve_fit"
        step_speed = 10000 #making the erf-function a almost step-function
        return height + size * erf(step_speed*(time-center))
    best_fit = curve_fit(erf_background,time_data,func_data)
    return best_fit

#constrained exponential function
def expfunc(time, size, height, slope):
    time = np.array(time)
    return height+np.exp(-slope*(time-size))-np.exp(0)

#curvefit constrained exp function
def fit_exp_background(time_data, func_data, center):
    time_data=np.array(time_data)
    func_data=np.array(func_data)
    best_fit = curve_fit(expfunc,time_data,func_data,p0=[center,func_data[0],0.005],bounds=([time_data[0], 0, 0],[time_data[len(time_data)-1],5000.,2.]))
    return best_fit

#save data as CSV
def save_to_csv(data, output_path=None, data_filename=None):
    """Converts data to dataframe (if necessary) and exports  as csv (comma sep)"""
    if type(data) is not pd.DataFrame:
        df = pd.DataFrame(data)
    else:
        df = data
    if data_filename == None:
        data_filename = input("Enter a name for the datafile")
    if output_path == None:
        df.to_csv(data_filename + '.csv', na_rep='NULL')
    else:
        df.to_csv(output_path + '\\' + data_filename + '.csv', na_rep='NULL')
    print("Data saved as " + data_filename)

#fit background and integrate spectrum
def integrate_spectrum(sequence,injection,detector,info): #detectors are "FID1A.ch" and "TCD2B.ch"
    spectrum = extract_data(sequence,injection,detector)
    results_array = {}
    error_array = {}
    
    for peak in info[detector]:
        mask = np.logical_and(spectrum[0]>info[detector][peak]['start'][0],spectrum[0]<info[detector][peak]['end'][0])
        peak_time = spectrum[0][mask]
        peak_count = spectrum[1][mask]
        background_mask_1 = np.logical_and(peak_time>info[detector][peak]['start'][0],peak_time<info[detector][peak]['start'][0]+info['settings'][detector]['background_range'][0])
        background_mask_2 = np.logical_and(peak_time>info[detector][peak]['end'][0]-info['settings'][detector]['background_range'][0],peak_time<info[detector][peak]['end'][0])
        peak_background_time = np.array([peak_time[background_mask_1],peak_time[background_mask_2]]).flatten()
        peak_background_count = np.array([peak_count[background_mask_1],peak_count[background_mask_2]]).flatten()
        
        erf_center = [i for i, j in enumerate(peak_count) if j == np.amax(peak_count)] #the step of the error-function is placed at peak-maximum
        erf_center = peak_time[erf_center[0]] #value instead of a vector
        
        background_fit = fit_erf_background(peak_background_time, peak_background_count, erf_center)
        #if background_fit[0][0] < -0.01: #and detector == FID_key: #negative steps will be given an exponential background
        #    background_fit = fit_exp_background(peak_background_time, peak_background_count, erf_center)
        #    renorm_peak_count = np.array(peak_count)-expfunc(peak_time,*background_fit[0])     
        #    area = (renorm_peak_count*np.mean(np.diff(peak_time))*60).sum() #60 is for seconds            
        #else: #If the exponential-background is not wanted, rewrite with only "else"-part of the loop
        renorm_peak_count = np.array(peak_count)-errorfunc(peak_time,erf_center,*background_fit[0])     
        area = (renorm_peak_count*np.mean(np.diff(peak_time))*60).sum() #60 is for seconds
        area_error = (peak_time[-1]-erf_center)*background_fit[0][0]
        
        results_array[peak] = abs(area)
        error_array[peak] = abs(area_error)
        
    for peak in info[detector]:
        for motherpeak in info[detector][peak]['mother_peak']:
            results_array[motherpeak] = results_array[motherpeak]-results_array[peak]
    return results_array, error_array


def integrate_spectrum_linear(sequence, injection, detector, info, name, save_plots=False, output_path = None):  # detectors are "FID1A.ch" and "TCD2B.ch"
    spectrum = extract_data(sequence, injection, detector)
    results_array = {}
    error_array = {}

    fig = plt.subplot()


    for peak in info[detector]:

        mask = np.logical_and(spectrum[0] > info[detector][peak]['start'][0],
                              spectrum[0] < info[detector][peak]['end'][0])
        peak_time = spectrum[0][mask]
        peak_count = spectrum[1][mask]
        #NOW instead of fitting an error function I want to fit a linear BG. even though not optimal, I think best is
        # a valley-to-valley approach with a separate consideration of peaks that lie on top of other peaks
        # this is critical for CH3OH&MeCHO, Acetone,EtCHO and acrolein -> here it would be optimal to fit an exponential
        #curve underneath if that looks good
        #algorithm borrowed from EC-MS analysis
        #here each peak is integrated separately in the for loops, not super efficient but works
        # finds the maximum in the defined region and then chooses start and end based on minima left and right of this maximum.
        # doesnt work if the peak is not significantly higher than the BG -> while loop added to take care of this
        data_slice = (peak_time, peak_count)
        # print(data_slice)
        start_slice = len(spectrum[1][spectrum[0] < info[detector][peak]['start'][0]])  # lines cut off when determining the index of the peakmax
        # start_slice = 0
        # print(start_slice)
        # print(tend)
        local_idx_peakmax = np.argmax(data_slice[1])
        # print("peakmax of " + str(data_slice[1][local_idx_peakmax]) + " at local index: " + str(local_idx_peakmax))
        # print(peak)
        print("Peakmax for " + str(peak) + " found at: " + str(data_slice[0][local_idx_peakmax]) + " min.")

        try:
            local_idx_min_start = np.argmin(data_slice[1][:local_idx_peakmax])

        except ValueError:
            print("Maximum value is at start or end of peak region for " + str(peak)+ " - double check peak limits!")
            print("Now trying to find peak maxium iteratively.")
            # print("initial length data slice=" + str(len(data_slice[1])))
            while local_idx_peakmax < 10:
                data_slice =(data_slice[0][int(len(data_slice[1])*0.01):-1],
                            data_slice[1][int(len(data_slice[1])*0.01):-1])
                # print(data_slice)
                # data_slice = new_slice
                # print("new length data slice=" + str(len(data_slice[1])))
                try:
                    local_idx_peakmax = np.argmax(data_slice[1])
                except ValueError:
                    print("No peakmax found")
                    break
                # print("local peak max index = " + str(local_idx_peakmax))
                # print(local_idx_peakmax < 10)
                # print(data_slice[0][local_idx_peakmax], data_slice[1][local_idx_peakmax])
            try:
                print("Peakmax found at: " + str(data_slice[0][local_idx_peakmax]) + " min.")
            except IndexError:
                data_slice = (peak_time, peak_count)
                local_idx_peakmax = 1
            local_idx_min_start = np.argmin(data_slice[1][:local_idx_peakmax])

        local_idx_min_end = np.argmin(data_slice[1][local_idx_peakmax:]) + local_idx_peakmax
        # print("local indizes" + str(local_idx_min_start) + " and " + str(local_idx_min_end))

        def_BG_idxlist = [0] + sorted([local_idx_min_start] + [local_idx_min_end]) + [-1]
        # print("BG indexlist: " + str(def_BG_idxlist))
        BG = np.interp(peak_time, peak_time[def_BG_idxlist], peak_count[def_BG_idxlist])

        #PLOT CURRENT PEAK plus background
        fig.plot(peak_time, peak_count, label=info[detector][peak]['name'], color = info[detector][peak]['color'] )
        fig.plot(peak_time, BG, ls=':', color = info[detector][peak]['color'])
        # fig.plot(peak_time, BG, label=info[detector][peak]['name'] + "BG", ls=':')

        #now for the actual integration part
        integral_peak = np.trapz(peak_count[local_idx_min_start:local_idx_min_end], peak_time[local_idx_min_start:local_idx_min_end]) *60
            #*60 to convert to seconds
        integral_peak_BG = np.trapz(BG[local_idx_min_start:local_idx_min_end], peak_time[local_idx_min_start:local_idx_min_end]) *60
        integral_BG_subtr = integral_peak - integral_peak_BG
        # print('For peak starting at t=' + str(spectrum[0][local_idx_min_start + start_slice]) + ': integral BG subtracted = ' + str(
        #     integral_BG_subtr) + " pA*s")

        #for gas injections: correct for the pressure as given in the data info or name from GC file

        try:
            pressure = sequence.injections[injection].metadata['pressure']
            print(pressure)
            area = integral_BG_subtr/pressure*1015
            print("Peak area normalized to atmospheric pressure.")
        except KeyError:
            area=integral_BG_subtr
            print("No pressure for correction found.")
        area_error=0

        results_array[peak] = abs(area)
        error_array[peak] = abs(area_error)

    #more details for plot
    time_frame = fig.get_xlim()
    mask =  np.logical_and(spectrum[0] > time_frame[0], spectrum[0] < time_frame[1])
    # mask =  np.logical_and(spectrum[0] > time_frame[0] + time_frame[0]*0.1, spectrum[0] < time_frame[1])
    x_spec,y_spec = spectrum[0][mask], spectrum[1][mask]
    fig.plot(x_spec, y_spec, "k:")
    fig.legend(loc='upper right')
    fig.set_xlabel("retention time / min")
    fig.set_ylabel("Signal intensity / pA")
    fig.axes.set_title(name)
    if save_plots == True:
        save_a_plot(output_path, '\\' + detector[0:4] +'_integrated_' + str(name), pdf = False)
    else:
        plt.show()


    for peak in info[detector]:
        for motherpeak in info[detector][peak]['mother_peak']:
            results_array[motherpeak] = results_array[motherpeak] - results_array[peak]
    return results_array, error_array