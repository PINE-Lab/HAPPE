# neuroscanio
Import .CNT, .EEG Neuroscan binary files as well as Neuroscan epoch file (.DAT) and Neuroscan event files (.EV2) into EEGLAB. There are also functions to import Neuroscan ASCII (text) location files from the command line as well as a beta function to export continuous CNT files from the command line.	

Note: In our experience, importing Neuroscan files is not an easy matter and no universal read function may exist. We include some other old tools below and hope they might help. Most of the documentation below is from 2002 but we included it as it might still be relevant for some researchers.

A collection of tools are included in [this data file](https://sccn.ucsd.edu/~arno/cntload.zip). They are detailed below.

The only function that I know which worked on my data (acquired with NEUROSCAN 3) were these 2 programs under WINDOWS/DOS-. 
However we do not have the source of these program and I was trying to read my data under UNIX.

CNTTOASC from this local folder (not available from Neuroscan anymore).
CNT2BIN from the neuroscan web site

OTHER PROGRAMS THAT READ CNT

CNT2ASC
(I put the link because the sources are no more accessible). However the sources do not work (output and empty file and say error in file format). I compiled them on Unix and Dos and it is the same.
but they should be close to the sources of  CNTTOASC. I tried to look at how they read the raw data. They use the channeloffset field to read blocks of data. However this field is 0 on my data.

LOADCNT of the neurotool software.

LDCNTB; Matlab function of Andrew James. Read well all dataset info but problem for the data itself (which are not continuous). This function is available on the previous page or down on this one.

ANALYSE_NEUROSCAN. Matlab function by J-R Duann, a co-worker of mine who have been working for Neuroscan (this guy knew what he was doing for sure). I did not test the function but it seams to me that it is only compatible with version. The function is here.

GNUROSCAN 
The header file gnuro.h which is included alongside cntcat.c in the
Gnuroscan distribution defines the offsets for relevant elements of
the Setup (S_), Electloc (EL_), and tagged-eeg (TEEG_) data structures.
Note that Neuroscan's .CNT format actually encompasses two formats, a
block-multiplexed one for data acquired from Neuroscan's own SYNAMPS
hardware, and a fully multiplexed one for non-SYNAMPS data. Gnuroscan
handles only the fully multiplexed (non-SYNAMPS) data, because I didn't
have access to any SYNAMPS data when I wrote the code. The data that
you mention in the final paragraph on this page, in blocks of forty
short integers, must be from SYNAMPS (info courtesy of Matthew Belmonte).

STAN SOFTWARE
no sources

DUMPCNT
cannot acess page

TESTREAD is small C program determining the offset of fields in the header. However, there is a problem with the size of the header (sethead.h) because offsets appear to be 4 bytes shifted (I redownloaded the sethead.h file from the Neuroscan site, but it's the same. I manually corrected for that but send me an email at arnosalk.edu (fake 'at' sign to prevent spam) if you figure out why.

eeg_load_scan_cnt from the EEG TOOLBOX (http://eeg.sourceforge.net). I couldn't get the function to work on my data.

After extensive searches, we realized that for our continuous CNT files, data was stored in blocks of 40 unsigned short integers for each channel [[40 bytes] * nb_channel] * nb_blocks. We couldn't find where this parameter (size of the block) was specified in the header, so we added the option 'blockread' and input the number 40 manually, ldcnt (then what I read is undistinguishable from the ASCII files I would read in Matlab after conversion by CNTTOASC or CNTTOBIN).

Time interval input for loadcnt function: If the imported data don't look like continuous EEG, try changing this number. Most often it should be 1 or 40, but other values may work. 

# Version history

Version 1.6
- Deleting readneurolocs to avoid conflicts

Version 1.5
- Fix issue for Octave compatibility

Version 1.4 update
- fix readneurolocs as it was importing an additional empty column

Version 1.3 update
- loadcnt can now read events for files larger than 1Gb. Contribution from edauer1 on Github.

Version 1.2 update
- pop_writeeeg.m was removed and placed back in the main EEGLAB distribution (it was a mistake to include it in the first place as it is designed to export EDF and BDF files, not Neuroscan files)
