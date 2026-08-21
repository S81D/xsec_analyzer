Guide for ANNIE Usage of STV-Analysis

--------------
***Standard Workflow***
--------------
1. Enter the ANNIE ToolAnalysis container: `./start_singularity.sh`
2. Set up the environment: `source setup.sh`
3. Prepare input files (if necessary): `root -l stvPrep.C`
4. Clean and compile: `make clean` and `make`
5. Univmake: `./univmake files_to_process.txt tutorial_bin_config.txt output.root`
6. XSec & Systematic Uncertainties: `./NC_analyzer`

-------------
***stvPrep***
-------------
This script (`stvPrep.C`) prepares ANNIE MC PhaseIITree ntuples for univmake. It creates and adds flag branches used in `tutorial_bin_config.txt` from existing branches to the ntuple file.
To use:
1. `hadd [target file] [source file 1] ...`
2. copy/move target file outside of `/pnfs/`
3. Edit target file path in `stvPrep.C` to match intended input
4. `root -l stvPrep.C`

**Notes**:
 - Input file path is hardcoded in `stvPrep.C` line 32
 - Input file gets updated, no existing branches are rewritten
 - Bug: DOES NOT RUN ON FILES LOCATED IN `/pnfs/` (files MUST be moved outside of `/pnfs/`)
 - DO NOT RUN DIRECTLY ON ANNIE MC PRODUCTION FILES (files get updated, make a copy first)
 - For the NC analysis, we use this to extract the XS + flux + MC stats uncertainties. Though it can handle data throughput, we have a separate XS extraction script (1D, flux-averaged).
 - Feel free to add extra flag/category branches! Please document any additions
 - Does not check if additional branches already exist
 - This version is different than the CC analysis workflow

------------------------------------------------
***Notable Differences between ANNIE & uBooNE***
------------------------------------------------
 - TTree Name - SystematicsCalculator.hh line 854 & univmake.C line 28 - TTree name for ANNIE files is 'phaseIITriggerTree'
 - POT Problem - SystematicsCalculator.hh line 861 - POT is hardcoded
  ANNIE MC ntuples do not include propagate the POT from GENIE through ToolAnalysis. POT is calculated externally for original GENIE files using POT files found here: `/pnfs/annie/persistent/users/doran/GENIE-016-empirical-model/genie/POT/`. Ensure POT on line 861 matches POT listed for onBNB in file_properties.txt. POT per event currently does not exist with ANNIE MC.
 - Spline weight - UniverseMaker.hh lines 89, 102 - Spline Weight is non-existent and unnecessary for ANNIE so the weight was replaced with 1.
 - Weight Types - UniverseMaker.hh lines 93-97 - ANNIE's flux weights are saved to individual branches, not to a "weight_flux_all"
 - ANNIE FV & flux - includes/AnnieGeometryTools.* - uses ANNIE FV and integrated flux window.

------------------------------------------------
***NCQE-specific workflow***
------------------------------------------------
I have archived this workflow to ensure the NCQE results are repeatable. This workflow differs from the CC / Joint analysis.

1. Assuming you have produced MC samples with weights, migrate the files (if needed) from `/pnfs/` to somewhere in `/exp/annie/data/` (ensuring to use `ifdh cp` if transferring from dCache). This is needed as we will directly append the TTrees. `/exp/annie/data/` is preferred as a storage location over `/exp/annie/app/` given the large volume of files.

2. Run `stvPrep.C` per the initial instructions under `***Standard Workflow***`. After execution, you will be left with N MC files that have been appended event selection and categories, with XS and flux weights appended together under an `all` branch.

4. Combine the MC files using `hadd` per the instructions under `***Standard Workflow***`. NOTE: as the world MC files are LARGE you will likely run into the 50 GB hadd limit if reproducing this workflow. The solution is to hadd the MC files in batches (assuming there are ~4000 files for the WORLD samples):
  - 'hadd [target file] MC_1*.root'
  - 'hadd [target file] MC_2*.root'
  - 'hadd [target file] MC_3*.root'
  - 'hadd [target file] MC_4*.root MC_5*.root MC_6*.root MC_7*.root MC_8*.root MC_9*.root MC_0*.root'

5. (assuming you hadd the files individually as described above) Because we see less dirt events in the data as compared to the simulation, we can "downsample" the dirt events to strip them from the MC record. We do this for each hadd file via: `root -l downsample_dirt.C` (edit the paths accordingly; currently the dirt scaling factor is set to alpha = 0.27). This script will omit (1-alpha) dirt events, and set all cross section weights for the dirt events = 1.0 (since they were scaled in the sideband). This step can be omitted or changed depending on the future choices of analyzers. James' CC-analysis handling of dirt events is done differently.

6. With the files downsampled (and now smaller than before), we can hadd them together to form one MC file with the appropriate weights: `hadd [target file] downsampled_MC_files*.root`

7. For FakeData / MC cross section extraction, `xsec_analyzer` requires one MC weighted file, one "Fake Data" MC file, and one "off-beam" file. For both the FakeData and offbeam file (assuming you are just running GENIE samples), we can run: `root -l keep_one_entry.C` which produces a single file with a single MC entry --> this can be used as the "off-beam" MC sample. We have a separate analysis workflow for determining the true off-beam sample for the actual cross section extraction in data. Update `file_properties.txt` and `files_to_process.txt` accordingly with the off-beam file. Set the triggers to 1 (since its just a single event).

8. We then can run `root -l stripTree.C` to strip out a majority of the branches to prepare for univmake. This should be done over the off-beam MC and to produce a separate FakeData file (that is a stripped copy of the full MC weighted file). The script will importantly leave in the CV tuning weights and thats necessary for the systematic calculators.

9. At this point you should now be left with 3 files:
 - off-beam/bkg file (say MC_bkg.root) that contains a single entry, CV weight, and NCQE event selection and truth branches (all other branches stripped)
 - MC fully weighted file, containing all MC events hadd-d together from the MC individual files that were processed through `stvPrep.C`
 - MC "Fake-Data" file that is the file above with most of the branches stripped, with only the CV weights and NCQE event selection and truth branches (all other branches stripped)

10. Update the paths in `files_to_process.txt`:
```
/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/MC_bkg.root                                                          # bkg file
/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled.root                # MC weighted file
/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled_stripped.root       # MC Fake-Data
```

11. Update the paths, POT, and trigger information in `file_properties.txt`:
```
# Paths to relevant files; column information: path_to_file, run_number(trivial for annie, set to 1), tag_name, triggers, POT
# Beam-off & MC should be assigned 0 POT, all runs should be set to 1
# As the Beam-off is MC, just assign 1 trigger
#
# Beam-off data (fake, MC single entry file)
/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/MC_bkg.root 1 extBNB 1 0
#
# Full beam-on data (MC FakeData)
/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled_stripped.root 1 onBNB 2358240 4.417e20
#
# MC samples
/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled.root 1 numuMC
```
The total trigger counts can be obtained by counting the total entries in the TTree. A macro has been provided to extract the total counts (triggers), as well as the number of signal/background events (true and reconstructed) + the efficiency: `CheckEventRates.C`, executed via: `root -l CheckEventRates.C`. 

12. After ensuring the correct MC POT is hardcoded in SystematicsCalculator.hh and that your binning is set in `tutorial_bin_config.txt`, run univmake: `./univmake files_to_process.txt tutorial_bin_config.txt output.root`. After running (may take a while) you will be left with `output.root`.

13. For the full breakdown of the systematics, run `./NC_analyzer`.

---------------------------
***Most up to date files***
---------------------------
 - BG file - (just a regular ntuple with only a single event)
`/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/MC_bkg.root`

 - MC files - File name, ~# of events, ~POT
`/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled.root 2358240 4.417e20`

 - fake data - Ntuples with only CV weight
`/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled_stripped.root 2358240 4.417e20`

 - stv output - Currently ran with MC and fake beam data
`/exp/annie/data/users/doran/xsec_analyzer_files/output.root`

