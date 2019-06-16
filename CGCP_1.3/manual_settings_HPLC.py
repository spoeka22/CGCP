"""
This script holds the information needed in 'CGCP' in order to treat data.
"""

"""
MANUEL SETTINGS SECTION
"""
#Path to files. The path here, can be both relative and absolute
path = r"C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\HPLC product analysis\ANNA001 2019-05-23 13-45-07"

output_path = r'C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\HPLC product analysis\test'

save_plots = False

#settings for calculation of concentration and electrochemical characteristics

calibration_path = r"C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\Calibration\new GC calibration\calibration_lines_HSGC150219+GC200519.csv"
echem = False #calculate FE and partial current density
exp_length = 30 #CA duration in min
echem_info = {"standard": {'liquid_volume': 12, 'total_charge': 0.1, 'electrode_area': 2},
             "AW Au 005": {'liquid_volume': 12.85, 'total_charge': 0.033512166, 'electrode_area':16.9}, #should be 006, but for gas products wrong name
"AW Au 007": {'liquid_volume': 12.65, 'total_charge': 0.114965769099999, 'electrode_area': 15.9},
"AW Au 008": {'liquid_volume': 12.45, 'total_charge': 0.3871161945, 'electrode_area': 17.9},
"AW AuPd 9010 004": {'liquid_volume': 12.35, 'total_charge': 0.5205460207, 'electrode_area': 15.9},
"AW AuPd 9010 005": {'liquid_volume': 12.77, 'total_charge': 0.004187647997, 'electrode_area': 16.8},
"AW AuPd 9010 006": {'liquid_volume': 12.7, 'total_charge': 0.5023260762534, 'electrode_area': 15.5},
"AW AuPd 1090 003": {'liquid_volume': 12.5, 'total_charge': 0.4226389751, 'electrode_area': 55},
"AW AuPd 9010 007": {'liquid_volume': 12.7, 'total_charge': 0.0430297937, 'electrode_area': 37.3},
"AW AuPd 1090 004": {'liquid_volume': 12.9, 'total_charge': 0.06137843503, 'electrode_area': 16.8},
"AW AuPd 1090 005": {'liquid_volume': 12.65, 'total_charge': 0.0438387044, 'electrode_area': 88.0},
"AW AuPd 1090 006": {'liquid_volume': 12.55, 'total_charge': 0.226526717, 'electrode_area': 38.4},

              }#total electrolyte volume in mL, total charge at the end of the CA experiment in C, electrode area cm2
electrons = {} #dict to be filled below for electron transferred for each compound.


#keys for detectors, which will depend on the actual GC
FID_key = 'rid1A.ch'
# TCD_key = None
TCD_key = 'dad1E.ch'

detector_keys = [FID_key, TCD_key]

#deactivated in script
# #injections to be skipped at the beginning
# skip_start = 0 #set to 1, and it will skip the first injection. The injection-times are unaltered
#
# #injections to be skipped at the end
# skip_end = 0 #set to 1, and it will skip the last injectionplt.show

#A lot of plot-settings
import matplotlib as mpl

#Titles written in figures
plots_title = ''

#plot-design
font = {'family' : 'serif',#'Palatino Linotype',
        #'weight' : 'bold',
        'size'   : 18}
mpl.rc('font', **font)

#Figure size
from pylab import rcParams
rcParams['figure.figsize'] = 10, 7
            
"""
INFO FOR BACKGROUND FITTING AND INTEGRATION

Names ('keys') and ordering of the peaks within the fit_info is abitrary, since fit_info is a dictionary.
The disignation 'P1', 'P2' and so on are not used for anything. The gc_parser code will just iterate over
all entries in the dictionary anyway.

New peaks can be added by simply 'copy-paste' a section, containing a peak. It is important to have different
entries (P1,P2 etc.), else earlier entries with same name will be overwritten when integrating spectra.

Peaks that need to be excluded, can be out-commented by # or sectional wise, like 'P2' in TCD.

Compound-names are written in a latex-syntax.

Colors are used in the plots
"""
#for the test_data
fit_info = {}
fit_info['settings'] = {}
fit_info['settings'][FID_key] = {}
fit_info['settings'][FID_key]['background_range'] = [0.05] #minutes
# fit_info['settings'][TCD_key] = {}
# fit_info['settings'][TCD_key]['background_range'] = [0.05] #minutes

fit_info[FID_key] = {}

#gas products:
#FID: CO, methane, CO2, ethylene
#TCD: H2 (4.2), air (N2+O2) (5.x)



fit_info[FID_key]['CO'] = {}
fit_info[FID_key]['CO']['start'] = [1.25]
fit_info[FID_key]['CO']['end'] =  [1.35]#[4]
fit_info[FID_key]['CO']['name'] = 'CO'
fit_info[FID_key]['CO']['color'] = 'grey'
fit_info[FID_key]['CO']['mother_peak'] = []
electrons['CO'] = 4

#
fit_info[FID_key]['CH4'] = {}
fit_info[FID_key]['CH4']['start'] = [1.35]
fit_info[FID_key]['CH4']['end'] =  [1.5]#[4]
fit_info[FID_key]['CH4']['name'] = 'CH4'
fit_info[FID_key]['CH4']['color'] = 'r'
fit_info[FID_key]['CH4']['mother_peak'] = ['CO']
electrons['CH4'] = 0 #not a product?!

#
fit_info[FID_key]['CO2'] = {}
fit_info[FID_key]['CO2']['start'] = [1.72]
fit_info[FID_key]['CO2']['end'] =  [2.1]#[4] #2.15 is A BIT TOO FAR FOR for propene oxidation results
fit_info[FID_key]['CO2']['name'] = 'CO2'
fit_info[FID_key]['CO2']['color'] = 'brown'
fit_info[FID_key]['CO2']['mother_peak'] = []
electrons['CO2'] = 6
#
fit_info[FID_key]['CH2CH2'] = {}
fit_info[FID_key]['CH2CH2']['start'] = [2.4]
fit_info[FID_key]['CH2CH2']['end'] =  [2.6]#[4]
fit_info[FID_key]['CH2CH2']['name'] = 'CH2CH2'
fit_info[FID_key]['CH2CH2']['color'] = 'g'
fit_info[FID_key]['CH2CH2']['mother_peak'] = []
electrons['CH2CH2'] = 0 #not a product?!
#

# fit_info[FID_key]['CH2CH2CH3'] = {}
# fit_info[FID_key]['CH2CH2CH3']['start'] = [2.10]
# fit_info[FID_key]['CH2CH2CH3']['end'] =  [2.4]#[4]
# fit_info[FID_key]['CH2CH2CH3']['name'] = 'propene?'
# fit_info[FID_key]['CH2CH2CH3']['color'] = 'grey'
# fit_info[FID_key]['CH2CH2CH3']['mother_peak'] = []
# electrons['CH2CH2CH3'] = 0 #not a product?!
#
# #what if this is formaldehyde?? then it would be 2e- per formaldehyde, 6e- per propene
#
#
# fit_info[FID_key]['MeOH'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[FID_key]['MeOH']['start'] = [2.9]
# fit_info[FID_key]['MeOH']['end'] = [3.24]
# fit_info[FID_key]['MeOH']['name'] = 'MeOH'
# fit_info[FID_key]['MeOH']['color'] = 'orange'
# fit_info[FID_key]['MeOH']['mother_peak'] = []
# electrons['MeOH'] = 0 #formally no oxidation/reduction per whole propene molecule
# #
#
# fit_info[FID_key]['MeCHO'] = {}
# fit_info[FID_key]['MeCHO']['start'] = [3.20]
# fit_info[FID_key]['MeCHO']['end'] =  [3.6]#[4]
# fit_info[FID_key]['MeCHO']['name'] = 'MeCHO'
# fit_info[FID_key]['MeCHO']['color'] = 'b'
# fit_info[FID_key]['MeCHO']['mother_peak'] = []
# electrons['MeCHO'] =  2 #(correct??)
# #
# fit_info[FID_key]['EtOH'] = {}
# fit_info[FID_key]['EtOH']['start'] = [5.6]
# fit_info[FID_key]['EtOH']['end'] = [6.2] #6.5 for standards
# fit_info[FID_key]['EtOH']['name'] = 'EtOH'
# fit_info[FID_key]['EtOH']['color'] = 'brown'
# fit_info[FID_key]['EtOH']['mother_peak'] = []
# electrons['EtOH'] = 0 # formally no oxidation (or reduction)
# #
# fit_info[FID_key]['Unknown01'] = {}
# fit_info[FID_key]['Unknown01']['start'] = [6.2]
# fit_info[FID_key]['Unknown01']['end'] = [6.7] #6.5 for standards
# fit_info[FID_key]['Unknown01']['name'] = 'Unknown01'
# fit_info[FID_key]['Unknown01']['color'] = '#7d40c6'
# fit_info[FID_key]['Unknown01']['mother_peak'] = []
# electrons['Unknown01'] = 0
# #
# fit_info[FID_key]['Acrolein'] = {}
# fit_info[FID_key]['Acrolein']['start'] = [6.7]
# fit_info[FID_key]['Acrolein']['end'] = [7.1] #7.1 for samples with low acrolein (Au, Au rich), 8.5 for high acrolein
# fit_info[FID_key]['Acrolein']['name'] = 'Acrolein?'
# fit_info[FID_key]['Acrolein']['color'] = '#40c69b'
# fit_info[FID_key]['Acrolein']['mother_peak'] = []
# electrons['Acrolein'] = 4
# #
# fit_info[FID_key]['EtCHO'] = {}
# fit_info[FID_key]['EtCHO']['start'] = [7.25] #[7.15] for low acrolein ([7.25] for high)
# fit_info[FID_key]['EtCHO']['end'] =  [7.8] # [9.1] for low acrolein ([7.65] for high)
# fit_info[FID_key]['EtCHO']['name'] = 'EtCHO'
# fit_info[FID_key]['EtCHO']['color'] = 'g'
# fit_info[FID_key]['EtCHO']['mother_peak'] = [] #["Acrolein"] #deactivate for low acrolein concentrations
# electrons['EtCHO'] = 2
# #
# fit_info[FID_key]['Acetone'] = {}
# fit_info[FID_key]['Acetone']['start'] = [7.85]
# fit_info[FID_key]['Acetone']['end'] = [9.0] #9.0 for higher acetone concentration 8.5 for low acetone
# fit_info[FID_key]['Acetone']['name'] = 'Acetone'
# fit_info[FID_key]['Acetone']['color'] = 'r'
# fit_info[FID_key]['Acetone']['mother_peak'] = [] #["Acrolein"] #deactivate for low acrolein concentrations
# electrons['Acetone'] = 2
# #
# fit_info[FID_key]['IPA'] = {}
# fit_info[FID_key]['IPA']['start'] = [10.6]
# fit_info[FID_key]['IPA']['end'] = [13]
# fit_info[FID_key]['IPA']['name'] = 'IPA'
# fit_info[FID_key]['IPA']['color'] = 'purple'
# fit_info[FID_key]['IPA']['mother_peak'] = []
# electrons['IPA'] = 0 #not echem
# #
# fit_info[FID_key]['Unknown02'] = {}
# fit_info[FID_key]['Unknown02']['start'] = [18.1]
# fit_info[FID_key]['Unknown02']['end'] = [19.6] #6.5 for standards
# fit_info[FID_key]['Unknown02']['name'] = 'Unknown02'
# fit_info[FID_key]['Unknown02']['color'] = 'm'
# fit_info[FID_key]['Unknown02']['mother_peak'] = []
# electrons['Unknown02'] = 0
# #
# fit_info[FID_key]['Unknown03'] = {}
# fit_info[FID_key]['Unknown03']['start'] = [40]
# fit_info[FID_key]['Unknown03']['end'] = [44] #6.5 for standards
# fit_info[FID_key]['Unknown03']['name'] = 'Unknown03'
# fit_info[FID_key]['Unknown03']['color'] = 'c'
# fit_info[FID_key]['Unknown03']['mother_peak'] = []
# electrons['Unknown03'] = 0
# #
# fit_info[FID_key]['Unknown04'] = {}
# fit_info[FID_key]['Unknown04']['start'] = [45]
# fit_info[FID_key]['Unknown04']['end'] = [50] #6.5 for standards
# fit_info[FID_key]['Unknown04']['name'] = 'Unknown04'
# fit_info[FID_key]['Unknown04']['color'] = '#780000'
# fit_info[FID_key]['Unknown04']['mother_peak'] = []
# electrons['Unknown04'] = 0


fit_info[TCD_key] = {}

fit_info[TCD_key]['H2'] = {}
fit_info[TCD_key]['H2']['start'] = [4.15]
fit_info[TCD_key]['H2']['end'] =  [4.55]#[4]
fit_info[TCD_key]['H2']['name'] = 'H2'
fit_info[TCD_key]['H2']['color'] = 'b'
fit_info[TCD_key]['H2']['mother_peak'] = []
electrons['H2'] = 2 #not a product?!


fit_info[TCD_key]['air'] = {}
fit_info[TCD_key]['air']['start'] = [5.18]
fit_info[TCD_key]['air']['end'] =  [5.8]#[4]
fit_info[TCD_key]['air']['name'] = 'air'
fit_info[TCD_key]['air']['color'] = 'grey'
fit_info[TCD_key]['air']['mother_peak'] = []
electrons['air'] = 0 #not a product?!


fit_info[TCD_key]['CO_TCD'] = {}
fit_info[TCD_key]['CO_TCD']['start'] = [8.07]
fit_info[TCD_key]['CO_TCD']['end'] =  [8.5]#[4]
fit_info[TCD_key]['CO_TCD']['name'] = 'CO_TCD'
fit_info[TCD_key]['CO_TCD']['color'] = '0.75'
fit_info[TCD_key]['CO_TCD']['mother_peak'] = []
electrons['CO_TCD'] = 0 #not a product?!

# #
#
# fit_info[TCD_key]['CO2'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['CO2']['start'] = [3.6]
# fit_info[TCD_key]['CO2']['end'] = [4.5]
# fit_info[TCD_key]['CO2']['name'] = 'CO$_2$/3.8 min'
# fit_info[TCD_key]['CO2']['color'] = 'b'
# fit_info[TCD_key]['CO2']['mother_peak'] = []
#