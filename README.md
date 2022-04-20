Copyright 2022 Dmitri Kosenkov

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


Installation Instructions

For Windows10:

	1.	Download and install 7-Zip https://www.7-zip.org/download.html

	2.	Download PyFREC: https://github.com/DKosenkov/PyFREC/archive/refs/heads/main.zip

	3.	Open and extract (unpack) PyFREC-main.zip with 7-Zip

	4.	Download and install Anaconda3 https://www.anaconda.com/

	5.	Start Anaconda Prompt (Anaconda3) form the start menu

	6.	Update conda:

	conda update conda
	7.	Create conda environment and install required packages for PyFREC :
	conda create --name pyfrec python=3.8 numpy=1.21.5 scipy=1.7.3 py=1.11.0 pytest=6.1.1

	8.	Activate environment for PyFREC:
	conda activate pyfrec

	9.	Change directory and run tests (default download directory: C:\Users\<User Name>\Downloads\). Then, inspect tests.log for potential warning and errors.
	cd C:\Users\dkosenko\Downloads\PyFREC-main\src\tests
	tests.bat

	10.	Change directory and run sample calculations. This will create “.log” files as output.
	Sample fluorescence lifetime calculations:
	cd C:\Users\dkosenko\Downloads\PyFREC-main\samples\00_fl_lifetimes
	run.bat

	Sample Forster calculations:
	cd C:\Users\dkosenko\Downloads\PyFREC-main\samples\01_bodipy_dte\Forster
	run.bat
	cd C:\Users\dkosenko\Downloads\PyFREC-main\samples\02_bodipy_trz\Forster
	run.bat

	11.	PyFREC usage in Anaconda3 prompt:

	cd <directory with your input files>
	conda activate pyfrec
	C:\Users\dkosenko\Downloads\PyFREC-main\pyfrec.py -f input_file.ini
	or
	C:\Users\dkosenko\Downloads\PyFREC-main\pyfrec.py -f  .\input_file.ini 1> .\output_file.log 2>&1

For Linux (Mac OS):

	1.	Download PyFREC:
	wget https://github.com/DKosenkov/PyFREC/archive/refs/heads/main.zip

	2.	Unpack PyFREC:
	unzip -a main.zip

	3.	Download and install Anaconda3 https://www.anaconda.com/
	
	4.	Update conda and initialize it (you may need to log-off and log-in again after initialization):
	conda update conda
	conda init bash
	
	5.	Create conda environment and install required packages for PyFREC:
	conda create --name pyfrec python=3.8 numpy=1.21.5 scipy=1.7.3 py=1.11.0 pytest=6.1.1
	
	6.	Activate environment for PyFREC:
	conda activate pyfrec

	7.	Change directory to the location of the unpacked PyFREC (e.g. /home/<User Name>):
	cd /home/dmytro/main

	8.	Change directory and run tests (e.g. download directory: /home/<User Name>/PyFREC-main). Then, inspect tests.log for potential warning and errors.
	cd /home/dmytro/PyFREC-main/src/tests
	chmod +x tests.sh
	./tests.sh

	9.	Change directory and run sample calculations. This will create “.log” files as output.
	Sample fluorescence lifetime calculations:
	cd /home/dmytro/PyFREC-main/samples/00_fl_lifetimes
	chmod +x run.sh
	run.sh
	
	Sample Forster calculations:
	cd /home/dmytro/PyFREC-main/samples/01_bodipy_dte/Forster
	chmod +x run.sh
	./run.sh

	cd /home/dmytro/PyFREC-main/samples/02_bodipy_trz/Forster
	chmod +x run.sh
	./run.sh

	10.	PyFREC Usage

	cd <directory with your input files>

	conda activate pyfrec
	/home/dmytro/PyFREC-main/pyfrec.py -f input_file.ini
	or
	/home/dmytro/PyFREC-main/pyfrec.py -f  ./input_file.ini 1> ./output_file.log 2>&1

	conda deactivate pyfrec
	
	Parameters: 
	input_file.ini  - main input file 
	output_file.log  - main output file
	
	(See user manual for additional details)
