     _____ _____ _____ _____ ____  _____ _____ _____ 
    |     |  _  |   __|   | |    \|   __|     |_   _|
    |  |  |   __|   __| | | |  |  |   __|   --| | |  
    |_____|__|  |_____|_|___|____/|_____|_____| |_|  
# openDECT_Light
A lighter version of OpenDECT
Here are described the major steps in the code:

1. Read CT scan DICOM file using the pydicom library
2. Select an area of interest, typically where rock samples were taken for core analysis.
3. Calculate density based on calibration.
4. Derive porosity, mean atomic number and photoelectric factor for each pixel in this region.
##Dependencies
- PyQt4
- [pydicom] (http://pydicom.readthedocs.io/en/stable/getting_started.html)
- Tcl/Tk (http://www.tcl.tk/software/tcltk/8.6.html)
