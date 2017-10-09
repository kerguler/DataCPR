# DataCPR
Continuous Plankton Recorder (CPR) Data Analyser

## Installation
The R package can be installed from the command line,

    R CMD install DataCPR_x.x.tar.gz

to be loaded easily at the R command prompt.

    library(DataCPR)

## Usage
In order to read the datasets into the dedicated class, we can use the following command:

    tow <- towClass(tow.log="data/9_TowLog_R.csv",
                    tow.ais="data/9_TowAIS_R.csv",
                    tow.ctd="data/9_CTD_R.csv",
                    tow.pci="data/9_PCI_R.csv",
                    tow.id=9,
                    silk.start=0,
                    silk.end=26.9,
                    time.offset=-2)

Then, we can view the summary of the class contents with the followings:

    print(tow)

    print(tow@data.log)
    print(tow@data.ais)
    print(tow@data.ctd)
    print(tow@data.pci)

    plot(tow)

Please note that, as of this moment, the package is under construction, and several functionalities are missing. For more information, please contact me at [k.erguler@cyi.ac.cy](mailto:k.erguler@cyi.ac.cy).
