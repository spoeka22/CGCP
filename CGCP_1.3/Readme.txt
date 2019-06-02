----------------------------------------------------------------------------------------------
   ______      ______      ______       __       _______      _______  _________   _______    
 .´ ___  |   .´ ___  |    |_   __ \    /  \     |_   __ \    /  ___  ||_   ___  | |_   __ \   
/ .´   \_|  / .´   \_|      | |__) |  / /\ \      | |__) |  |  (__ \_|  | |_  \_|   | |__) |  
| |    ____ | |             |  ___/  / ____ \     |  __ /    `.___`-.   |  _|  _    |  __ /   
\ `.___]  _|\ `.___.´\     _| |_   _/ /    \ \_  _| |  \ \_ |`\____) | _| |___/ |  _| |  \ \_ 
 `._____.´   `._____.´    |_____| |____|  |____||____| |___||_______.´|_________| |____| |___|
----------------------------------------------------------------------------------------------
For use at DTU Physics                                                    GC Parser 1.2 (2018)

CGCP - Continuous Gas Chromatograph Parser

(+) Initial notes
- This code has been written to work well, launched from a commando-prompt (so not spyder or 
any sort). Wether it works well on Spyder is up to the user, and (s)he is free to rewrite the 
code for this purpose.
- The logo is from before 1.1, were the very general name GC_parser was used.

(+) Requirements
- Working Python 3.6 distribution with the following packages: numpy, matplotlib, scipy, unicodecsv
- That your data is ordered in the right way (From agilient gas chromatographs). The folder 
  "test_data" is filled with "good" data.
- If you have the package PyExpLabSys from github, you should read the lines 64-66 in "CGCP_1.2.py".
  Else a good tip is to drop the whole "Resources" folder into your python-path.
- 

(+) Running
- 	Make sure that both "CGCP_1.2.py", "manual_settings.py" and a folder called "Resources" is
	at the same directory.

- 	Set the values in "manual_settings.py", which is a script that is exicuted from the main 
	script. Especially the paths for loading and saving data is important to take note on.

- 	Run "CGCP_1.2.py" with python from a commando-prompt (fx. Anaconda Prompt, or any prompt 
	that can launch python 3.0+). The syntax would be "python CGCP_1.2.py" or "python3 CGCP_1.2.py".

- 	use "plt.show()" to print figures. The matplotlib GUI can be used to zoom in on the figures.
	The two interactive displays contains sliding bars in the buttom for scrolling through the
        raw data. In the case something is off, change the specific limits in "manual_settings.py" 
	for the "fit_info" dictionary.

(+) Help and questions
Write to thoe@physics.dtu.dk or find me in 307:044 at DTU. 

(+) Thanks
To Kenneth Nielsen for the PyExpLabSys package available on github. The whole package is not used, but the script "chemstaion.py" is written for this package.

—————————————————



Version history of gc_parser.py



0.1.0

- Load temperature and gc-data
- full raw data rainbow plot
- FID and TCD-signal plot based on area integration from Agilent software



0.2.0 (Faulty)

- fitting using "skewed gaussian”. By curtesy of Simon.



0.2.1

- fitting using "exponential modelled gaussian”.



0.2.2

- added constraints of fit.

- added better start-guesses for fit.



0.3.0

- numerical integration.

- background fitted with erf-function (or erfc)

- multiple spectra could be treated in sequence.

- FID and TCD-signal plot based on numerical integration



0.4.0

- final written code

- functions and fit_info moved to separate script

- Added constrained exponential background, for when peaks in spectrum were placed in close vicinity. “Constrained” in the sense that it has a lesser degree of freedom but is more numerical stable.



1.0

- identical in many ways to version 0.4.0.

- included GC-calibrated values of FID and TCD detector for plotting real concentrations



1.1

- fixed time
- stamps for GC-data, so they are loaded from metadata instead of being generated

- changed temperature-dependant measurements, so point taken at temperature
- changes beyond 4 degrees are excluded

Version history of CGCP.py

1.2
- Major rewriting to shorten the code.
- fixed time for real, compared to 1.1
- changed loaded functions to actual imported functions
- moved often changed settings to a seperate script
- wrote the script, so PyExpLabSys functions are included in the "Resources"

1.3
- fixed bugs (still the time)
- included the interactive windows for scrolling through raw data
- Made the "light" version into the only one, accesible on CINF-WIKI.
