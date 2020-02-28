# Smith-Waterman-Implementation
Homework 1 for CBB 752

## Installation

### Installation using pip
~~~~~~~~~~~~~~~~
 pip install git+https://github.com/KevinBastianYang/SW.git
~~~~~~~~~~~~~~~~
### Installation using git
~~~~~~~~~~~~~~~~
git clone https://github.com/KevinBastianYang/SW.git
~~~~~~~~~~~~~~~~


## Usage

### Demostration: using sample files(under /files/ folder)
~~~~~~~~~~~~~~~~
from SW.SW import runSW

runSW(inputFile="sample-input.txt",scoreFile="blosum62.txt",openGap=-2,extGap=-1)
~~~~~~~~~~~~~~~~
### Output files will be under the location with name "SW_align_results.txt"
