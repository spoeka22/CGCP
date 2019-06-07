#This is script is specifically made for data-treatment of chromatography data sampled continously.
#
#First part will load data based on the package called PyExpLabSys, which can be found on github
#Thanks goes to the authors from CINF/PyExpLabSys (especially Kenneth Nielsen, that personally helped the author of this script).
#
#second part will integrate peaks in all spectra, estimate area and plot a time-series of
#the amounts from each peak. This fitting is based on the values given in the entity
#"fit_info", which holds general settings for the detectors (FID, TCD), and parameters
#for each peak for generating constraints. These parameters are written in manually based
#on the raw data from the first part of the script.
#
#Questions can be directed to Thomas Smitshuysen, thoe@physics.dtu.dk

#edited by Anna Winiwarter to make useful for automatic treatment of HS-GC data. OBS: important is the data format that
#the data is being saved in

import sys
from pprint import pprint  # pprint == pretty print

# Python version check
if sys.version_info.major < 3: #courtesy of Keni 
    print("No! Use Python 3")
    sys.exit(1)

"""
#START MANUAL SETTINGS
""" 
exec(open("manual_settings.py").read()) #This solution is generally 


"""
#END MANUAL SETTINGS
"""

print(' ')
print('All settings are configured:')
print('* Path of data: '+path)
print('* Keys for detectors in data-object:')
print('  FID: '+FID_key)
print('  TCD: '+TCD_key)
print(' ')

"""
#Load Functions
"""
import numpy as np
from matplotlib import pyplot as plt
import itertools
import pandas as pd

#from PyExpLabSys.file_parsers.chemstation import Sequence
#if PyExpLabSys in installed. Use line 64 instead of 66, and delete the files "cached_property.py","chemstation.py" and "supported_versions" within the folder "Resources"
from Resources.chemstation import Sequence

# from Resources.CGCP_utilities import extract_data
# from Resources.CGCP_utilities import integrate_spectrum
from Resources.CGCP_utilities import integrate_spectrum_linear
from Resources.CGCP_utilities import save_to_csv
from Resources.CGCP_utilities import extract_sample_info
from Resources.CGCP_utilities import save_a_plot
# from Resources.CGCP_utilities import all

# from Resources.upgrade_GC_data_to_series import create_csv_files

print(' ')
print('All functions are loaded.')
print(' ')

"""
#Load and plot raw GC data
"""

# Load Sequence. 
sequence = Sequence(path)

#make time-axis for plot of raw data (not retention time) => only used for making heatmap plot which doesnt make sense for me to use
# first_injection = sequence.injections[0].metadata['injection_date_timestruct']
# last_injection = sequence.injections[len(sequence.injections)-1].metadata['injection_date_timestruct']
# Time_start = first_injection[7]*24*60*60+first_injection[3]*60*60+first_injection[4]*60+first_injection[5]
# Time_end = last_injection[7]*24*60*60+last_injection[3]*60*60+last_injection[4]*60+last_injection[5]
# injection_spacing = (Time_end-Time_start)/len(sequence.injections) #average injection spacing, this might not be true


#instead either get concentration or experiment number from data
concentration_list, name_list, calibration, pressure_list = extract_sample_info(sequence)
#
# print(concentration_list)
# print(name_list)
# print(calibration)
# print(pressure_list)
# for injection in sequence.injections:
#     print(injection.metadata['pressure'])

def plot_vs_conc(data,title):
    fig, f = plt.subplots()
    f.plot(retention_time, data)
    f.set_xlabel('retention time [min]')
    f.set_ylabel('Intensity [a.u.]')
    f.set_title(title, loc='center', fontsize="small", wrap=True)
    # f.set_ylim(6,30)
    plt.tight_layout()
    # plt.draw()
    # plt.show()

#plotting of raw data for FID
if False:
#setup of raw data
    spectrum_length = len(sequence.injections[0].raw_files[FID_key].values)-10 #hardcoded to exclude the last 10 points, since the spectra-length vary this value will be the reference for the length, and generally the sampling is high, so a GC-spectrum could be 10000-100000 points long
    # print(spectrum_length)
    retention_time = sequence.injections[0].raw_files[FID_key].times[0:spectrum_length]
    # print(len(retention_time))
    if calibration == True:
        for injection, conc, number in itertools.zip_longest(sequence.injections, concentration_list, np.arange(len(concentration_list))):
            plot_vs_conc(injection.raw_files[FID_key].values[0:spectrum_length],'FID raw data ' + str(conc) + ' µM')
            if save_plots == True:
                save_a_plot(output_path, '\FID_raw_data_' + str(conc) + '_µM_inj_no_' + str(number), pdf=False)
    else:
        for injection, name in itertools.zip_longest(sequence.injections, name_list):
            # print('test_rawdata')
            plot_vs_conc(injection.raw_files[FID_key].values[0:spectrum_length], 'FID raw data ' + str(name))
            if save_plots == True:
                save_a_plot(output_path, '\FID_raw_data_' + str(name), pdf=False)

#now do the same for TCD data
if False:
    #setup of raw data
    spectrum_length = len(sequence.injections[0].raw_files[TCD_key].values)-10 #hardcoded to exclude the last 10 points, since the spectra-length vary this value will be the reference for the length, and generally the sampling is high, so a GC-spectrum could be 10000-100000 points long
    # print(spectrum_length)
    retention_time = sequence.injections[0].raw_files[TCD_key].times[0:spectrum_length]
    # print(len(retention_time))

    for injection, name in itertools.zip_longest(sequence.injections, name_list):
        plot_vs_conc(injection.raw_files[TCD_key].values[0:spectrum_length], 'TCD raw data ' + str(name))
        if save_plots == True:
            save_a_plot(output_path, '\TCD_raw_data_' + str(name), pdf= False)


print('Raw data has been plotted succesfully')
print(' ')

"""
#Data treatment
"""
if True:
#---comments below leftover from Thomas Smitshuysen
#Next part will setup the tools for fitting. haah just joking. Numerical integration for the win!
#though the background is fitted. But that fit could be done by anybody
#
#At some point, somebody should make a fit-function instead[1]. 
#[1] Kalambet et.al.: J. Chemometrics 2011; 25: 352–356
#
#The article is interesting, but the function itself is so numerically unstable that I wouldn't
#trust a gas-chromatography code, which fitted analytically. Ever! Unless I had written it ;-)

# functions for dealing with calibration lines:
    def calc_regression(x, y):
        """
        written by Lukas W. at some point...
        """
        coeff = np.polyfit(x, y, 1)
        p = np.poly1d(coeff)
        y_dach = p(x)
        abw = y - np.mean(y)
        abw2 = np.square(abw)
        b_gesch = y_dach
        abw_gesch = y - b_gesch
        abw_gesch2 = np.square(abw_gesch)
        R2 = 1 - (np.sum(abw_gesch2) / np.sum(abw2))
        print("k=%.3f, d=%.3f, R2=%.4f" % (coeff[0], coeff[1], R2))
        reg_lin_x = np.arange(0, max(concentration_list) + 50, 10)
        reg_lin_y = coeff[0] * reg_lin_x + coeff[1]
        return reg_lin_x, reg_lin_y, coeff[0], coeff[1], R2

    def plot_calibration_curve(result_dict, detector_key, title, ylabel, output_path, save_reg=True):
        fig, f = plt.subplots()
        reg_data = {"legend": ["k", "d", "R2"]}
        for entry in result_dict:
            print("Now plotting: " + str(entry))
            f.plot(concentration_list, result_dict[entry], ls="None", ms=10, marker="x",
                   label=fit_info[detector_key][entry]['name'], color=fit_info[detector_key][entry]['color'])
            # f.errorbar(injectiontimes_sec,result_dict[entry],yerr=error_dict[entry],fmt='o',capsize=1, color=fit_info[detector_key][entry]['color'])
            # but now we also want to get a calibration line if we're dealing with calibration data
            regression = calc_regression([0] + concentration_list,
                                         [0] + result_dict[entry])  # here the point 0/0 is added as a "fake" datapoint
            f.plot(regression[0], regression[1], ls=":", color=fit_info[detector_key][entry]['color'])
            if save_reg == True:
                reg_data[fit_info[detector_key][entry]['name']] = regression[2:]
        if save_reg == True:
            save_to_csv(reg_data, output_path)
        f.set_xlabel('concentration [µM]')
        f.set_ylabel(ylabel)
        f.tick_params('y')
        f.axes.set_title(title)
        f.legend(frameon=0)
        # f.set_title(plots_title,loc='right')
        plt.tight_layout()
        plt.draw()
        if save_plots == True:
            save_a_plot(output_path + '\calibration_curve', pdf=False)

    print('----------------------------------------------------')
    print(' ')
    print('Start integrating spectra.')

    integration_data = {}

    for detector in detector_keys:
        print('------------------------------------------------------')
        print("Now working on detector " + detector)
        print('------------------------------------------------------')
        # detector = detector[0:4]
        integration_data[detector + "raw"] = {}
        for peak in fit_info[detector]:
            integration_data[detector + "raw"][peak] = []
        if calibration == True:
            for injection_number, conc in itertools.zip_longest(range(0,len(sequence.injections),1), concentration_list):
                # print(f"{50*(injection_number-skip_start)/len(injectiontimes_sec):.1f} %\r",end="")
                name = str(conc) + " µM - injection no " + str(injection_number)
                spectrum_area = integrate_spectrum_linear(sequence,injection_number,detector, fit_info, name, save_plots, output_path)
                for peak in fit_info[detector]:
                    integration_data[detector + "raw"][peak].append(spectrum_area[0][peak])
                # plt.show()
        else:
            integration_data[detector + "raw"]['sample_info'] = name_list
            for injection_number, name in itertools.zip_longest(range(0,len(sequence.injections),1), name_list):
                spectrum_area = integrate_spectrum_linear(sequence, injection_number, detector, fit_info, name, save_plots, output_path)
                for peak in fit_info[detector]:
                    integration_data[detector + "raw"][peak].append(spectrum_area[0][peak])
                # plt.show()
            #save the integration data as csv
            save_to_csv(integration_data[detector + "raw"], output_path, name + detector + "_raw_integrals")

        #call to plot and calculate regression data (doesnt really make sense to have this as a separate file for HSGC data, but
        #might be useful if there's TCD data for gas analysis only
        if calibration == True:
            title = "GC " + str(detector)
            x_label = "Integrated signal [pA*s] " + str(detector)
            plot_calibration_curve(integration_data[detector + "raw"], detector, title, x_label, output_path)
        else:
            print("No calibration data plotted.")

    print('Integration of spectra is done')
    print(' ')

#read in a calibration dataset and calculate the concentration of the different substances
#also calculate the partial current densities and the faradaic efficiency based on input parameters
    for detector in detector_keys:
        # detector = detector[0:4]
        conc_eval = True
        try:
            calibration_path
        except NameError:
            print("No calibration path given. No calculations of concentration and echem possible.")
            conc_eval=False
        else:
            if calibration == False:
                integration_data[detector + "conc"] = {}
                reg_data = pd.read_csv(calibration_path)
                # print(reg_data)
                calibr_data = {}
                for item in list(reg_data.keys()):
                    calibr_data[item]= list(reg_data[item])
                del calibr_data['Unnamed: 0'] #not necessary but might prevent problems in iterations if I forget about it
                for peak in integration_data[detector + "raw"].keys():
                    if peak == "sample_info":
                        integration_data[detector + "conc"][peak] = integration_data[detector + "raw"][peak]
                    else:
                        try:
                            concentration = [(x - calibr_data[peak][1])/calibr_data[peak][0] for x in integration_data[detector + "raw"][peak]] #indices for k and d hardcoded now, it's not worth the effort to try and read that out of the list
                            integration_data[detector + "conc"][peak + "_conc"]=concentration
                        except KeyError:
                            print("No calibration data found for " + str(peak) + ", skipping concentration calculation.")
                        #need to write section where for each injection a bar-plot with concentrations of the individual
                        #compounts is made -> but this requires serious data restructuring -> maybe too much effort??
                save_to_csv(integration_data[detector + "conc"], output_path, integration_data[detector + "raw"]["sample_info"][-1] +  "_concentration")  # here it would be nice if it would actually also plot it

            # now use the electrochemistry data to calculate partial current density and FE
            if echem == True:
                echem_data = []
                for sequence_number in range(len(integration_data[detector + "conc"]['sample_info'])):
                    try:
                        echem_data.append(echem_info[integration_data[detector + "conc"]['sample_info'][sequence_number]])
                    except (NameError, KeyError):
                        echem_data.append(echem_info['standard'])
                # print("Echem data: ")
                # print(echem_data)
                integration_data[detector + "pcd"] = {}
                integration_data[detector + "fe"] = {}
                for peak in integration_data[detector + "raw"].keys():
                    if peak == "sample_info":
                        integration_data[detector + "pcd"][peak] = integration_data[detector + "raw"][peak]
                        integration_data[detector + "fe"][peak] = integration_data[detector + "raw"][peak]
                    else:
                        try:
                            # print("FID integration_data[detector + "conc"][peak] for peak: " + str(peak))
                            # print(integration_data[detector + "conc"][peak])
                            #input units: c: µmol/L, liquid volume: mL, exp length: min, electrode area: cm2  -> output: µA/cm2
                            pcd = [(c * 10 ** (-6) * ec_dat['liquid_volume'] * 10 ** (-3) * electrons[peak] * 96485 / (exp_length * 60 * ec_dat['electrode_area']) * 10 ** 6) for c, ec_dat in zip(integration_data[detector + "conc"][peak +"_conc"], echem_data)]
                            integration_data[detector + "pcd"][peak + "_pcd"] = pcd
                            fe = [c * 10 ** (-6) * ec_dat['liquid_volume'] * 10 ** (-3) * electrons[peak] * 96485 /  ec_dat['total_charge'] * 100 for  c, ec_dat in zip(integration_data[detector + "conc"][peak + "_conc"], echem_data)]
                            integration_data[detector + "fe"][peak + "_fe"]= fe
                        except KeyError:
                            print("No calibration data found for " + str(peak) + ", no echem data evaluated.")
                save_to_csv(integration_data[detector + "pcd"], output_path, integration_data[detector + "raw"]["sample_info"][-1] + "_partial_current_density")
                print("Partial current density saved: " + str(integration_data[detector + "pcd"]))
                save_to_csv(integration_data[detector + "fe"], output_path, integration_data[detector + "raw"]["sample_info"][-1] + "_FE")
                print("FE saved: " + str(integration_data[detector + "fe"]))

integration_data_to_save = {} #necessary to restructure dict to use the output function
for detector in detector_keys:
    integration_data_to_save.update(integration_data[detector + "raw"])
    if conc_eval == True:
      integration_data_to_save.update(integration_data[detector + "conc"])
      integration_data_to_save.update(integration_data[detector + "pcd"])
      integration_data_to_save.update(integration_data[detector + "fe"])
save_to_csv(integration_data_to_save, output_path, integration_data[detector + "raw"]["sample_info"][-1] + "_all")
                #need to write something to save this and plot this (and also in some way double check that what it does is actually correct - should be ok now)
"""
#END of script
"""

print('Program has finished')
print(' ')

# import code
# code.interact(local=locals())

# plt.show()