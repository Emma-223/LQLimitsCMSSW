# LQLimitsCMSSW
Scripts for computing limits, to be run inside a CMSSW area

## Installation instructions
First, make a CMSSSW area and install combine.
- See instructions: http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users

Next, `git clone` this package.

Finally, install the ruamel package (dependency from other scripts in rootNtupleMacrosV2 package):
`pip3 install --user ruamel.yaml`


## Usage instructions
### Nominal limit results
1. Estimate r-value scan range. This is done by running the AsymptoticLimits. Currently, this is done via batch jobs on condor.

Example: `python3 scripts/runLimits.py -d combCardFile.txt -n myLimitDir -e`

2. Read the estimation results into a json file.

Example: `python3 scripts/runLimits.py -d combCardFile.txt -n myLimitDir -e -r1

3. Use the estimatation results to submit limit-setting jobs.  Note that the r-value scan range estimations can be used from earlier results, as the scan ranges probably won't change much.

Example: `python3 scripts/runLimits.py -d combCardFile.txt -n myLimitDir -f myLimitDir/asymptoticLimits/condor/rValues.json`

4. Once the limit-setting jobs are done, extract the limits, make plots and tables:

Example: `python3 scripts/runLimits.py -d combCardFile.txt -n myLimitDir -f myLimitDir/asymptoticLimits/condor/rValues.json -r | tee myLimitDir/limitsNominal.log`


### Calculating limits in the beta vs. MLQ plane
Similar to the above, but specifying the `-b` flag to indicate that a scan over beta values should be done.
By default, this will also calculate the nominal limits.
One can include the flag `-s` to skip the nomimal limit calculation.
