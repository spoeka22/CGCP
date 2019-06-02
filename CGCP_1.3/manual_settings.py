"""
This script holds the information needed in 'CGCP' in order to treat data.
"""

"""
MANUEL SETTINGS SECTION
"""
#Path to files. The path here, can be both relative and absolute
# path = r"C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\GC product analysis\single runs"
# path = r"C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\Calibration\new GC calibration\data\ANNA_STANDARDS_NEWFLOW 2019-02-15 11-38-32 - Copy"

# path = r'C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\GC product analysis\ANNA_POR_NEW_FLOWRATE 2019-05-23 16-56-13'
# path = r'C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\GC product analysis\ANNA_POR_NEW_FLOWRATE 2019-05-24 18-35-58'
path = r'C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\GC product analysis\ANNA_POR_NEW_FLOWRATE 2019-05-27 17-56-38'
output_path = r"C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\201905_PdAufoams\GC product analysis\Data analysis"
save_plots = False

#settings for calculation of concentration and electrochemical characteristics
calibration_path = r"C:\Users\annawi\Desktop\Projects\Propene oxidation\Experiments\Calibration\new GC calibration\20190215_calibration\calibration_lines_HSGC150219.csv"
echem = True #calculate FE and partial current density
exp_length = 30 #CA duration in min
echem_info = {"standard": {'liquid_volume': 12, 'total_charge': 0.1, 'electrode_area': 2},
             "AW Au 006": {'liquid_volume': 12.85, 'total_charge': 0.033512166, 'electrode_area':16.9},
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

              }#total electrolyte volume in mL, total_charge = 0.1 #total charge at the end of the CA experiment in C, electrode area cm2
electrons = {} #dict to be filled below for electron transferred for each compound.


#keys for detectors, which will depend on the actual GC
FID_key = 'FID1A.ch'
# TCD_key = None
TCD_key = 'TCD2B.ch'

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

fit_info[FID_key]['MeOH'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
fit_info[FID_key]['MeOH']['start'] = [2.9]
fit_info[FID_key]['MeOH']['end'] = [3.24]
fit_info[FID_key]['MeOH']['name'] = 'MeOH'
fit_info[FID_key]['MeOH']['color'] = 'orange'
fit_info[FID_key]['MeOH']['mother_peak'] = []
electrons['MeOH'] = 0

fit_info[FID_key]['MeCHO'] = {}
fit_info[FID_key]['MeCHO']['start'] = [3.20]
fit_info[FID_key]['MeCHO']['end'] =  [3.6]#[4]
fit_info[FID_key]['MeCHO']['name'] = 'MeCHO'
fit_info[FID_key]['MeCHO']['color'] = 'b'
fit_info[FID_key]['MeCHO']['mother_peak'] = []
electrons['MeCHO'] = 0
#
fit_info[FID_key]['EtOH'] = {}
fit_info[FID_key]['EtOH']['start'] = [5.6]
fit_info[FID_key]['EtOH']['end'] = [6.5]
fit_info[FID_key]['EtOH']['name'] = 'EtOH'
fit_info[FID_key]['EtOH']['color'] = 'brown'
fit_info[FID_key]['EtOH']['mother_peak'] = []
electrons['EtOH'] = 0
#
fit_info[FID_key]['EtCHO'] = {}
fit_info[FID_key]['EtCHO']['start'] = [7.15]
fit_info[FID_key]['EtCHO']['end'] = [9.1]
fit_info[FID_key]['EtCHO']['name'] = 'EtCHO'
fit_info[FID_key]['EtCHO']['color'] = 'g'
fit_info[FID_key]['EtCHO']['mother_peak'] = []
electrons['EtCHO'] = 0
#
fit_info[FID_key]['Acetone'] = {}
fit_info[FID_key]['Acetone']['start'] = [7.88]
fit_info[FID_key]['Acetone']['end'] = [8.5]
fit_info[FID_key]['Acetone']['name'] = 'Acetone'
fit_info[FID_key]['Acetone']['color'] = 'r'
fit_info[FID_key]['Acetone']['mother_peak'] = ['EtCHO']
electrons['Acetone'] = 2
#
fit_info[FID_key]['IPA'] = {}
fit_info[FID_key]['IPA']['start'] = [10.6]
fit_info[FID_key]['IPA']['end'] = [13]
fit_info[FID_key]['IPA']['name'] = 'IPA'
fit_info[FID_key]['IPA']['color'] = 'purple'
fit_info[FID_key]['IPA']['mother_peak'] = ['EtCHO']
electrons['IPA'] = 0

fit_info[TCD_key] = {}
#
# fit_info[TCD_key]['CO2'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['CO2']['start'] = [3.6]
# fit_info[TCD_key]['CO2']['end'] = [4.5]
# fit_info[TCD_key]['CO2']['name'] = 'CO$_2$/3.8 min'
# fit_info[TCD_key]['CO2']['color'] = 'b'
# fit_info[TCD_key]['CO2']['mother_peak'] = []
#
# """
# fit_info[TCD_key]['P2'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['P2']['start'] = [5.6]
# fit_info[TCD_key]['P2']['end'] = [6]
# fit_info[TCD_key]['P2']['name'] = '5.8 min'
# fit_info[TCD_key]['P2']['color'] = 'gold'
# fit_info[TCD_key]['P2']['mother_peak'] = []
# """
#
# fit_info[TCD_key]['H2'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['H2']['start'] = [7.2]
# fit_info[TCD_key]['H2']['end'] = [7.65]
# fit_info[TCD_key]['H2']['name'] = 'H$_2$/7.4 min'
# fit_info[TCD_key]['H2']['color'] = 'c'
# fit_info[TCD_key]['H2']['mother_peak'] = []
#
# fit_info[TCD_key]['Ar'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['Ar']['start'] = [7.65]
# fit_info[TCD_key]['Ar']['end'] = [8.05]
# fit_info[TCD_key]['Ar']['name'] = 'Ar/7.8 min'
# fit_info[TCD_key]['Ar']['color'] = 'indigo'
# fit_info[TCD_key]['Ar']['mother_peak'] = []
#
# fit_info[TCD_key]['N2'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['N2']['start'] = [8.1]
# fit_info[TCD_key]['N2']['end'] = [8.7]
# fit_info[TCD_key]['N2']['name'] = 'N$_2$/8.2 min'
# fit_info[TCD_key]['N2']['color'] = 'm'
# fit_info[TCD_key]['N2']['mother_peak'] = []
#
# fit_info[TCD_key]['CO'] = {} #P is for peak in spectrum, and are ordered from low to high in retentiontime
# fit_info[TCD_key]['CO']['start'] = [9.1]
# fit_info[TCD_key]['CO']['end'] = [9.5]
# fit_info[TCD_key]['CO']['name'] = 'CO/9.2 min'
# fit_info[TCD_key]['CO']['color'] = 'saddlebrown'
# fit_info[TCD_key]['CO']['mother_peak'] = []
