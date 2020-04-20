# NanoAODv6 production instructions
Here follow indications on how to use these tools to produce standard nanoAOD version 6.


### Installation
Do the following in a CMSSW-supported environment.
```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src/
cmsenv
git cms-init
git remote add vrbouza https://github.com:/vrbouza/cmssw.git
git fetch vrbouza
git checkout CMSSW_10_2_18_modifiedLowPU
git cms-addpkg PhysicsTools/NanoAOD
scramv1 b -j 5
```


### Use
Everything happens inside the PhysicsTools/NanoAOD/test folder. Inside it, you have the multicrab.py file, with different options at the beginning that you should comment/uncomment. These options regard the datasets ready to be produced. The information for these MINIAOD datasets is stored in the datasets_MINIAOD folder in various json files: dataYEAR.json and mcYEAR.json store the dataset itself for each sample. Note that there is another classification level, up of the dataset: it is the dataset group, that allows you to easily join similar samples for similar processes together. In addition, in the nanoAODv6_cfg.json you can set the configurations for each one of the different sample cases that there are.

NOTE: for other kind of samples, although not yet tested, the possibility of using another "category" is available. E.g. one could use a susysignals.json configuration file in principle.

Once you have your dataset information configured, and the multicrab ready with the configuration, you can set the last CRAB details, if you want, in the CentralSettings.py file, inside the function LaunchCRABTask (here we refer to details such as DBS publication or not, storage site, and other CRAB details).

Afterwards, you can send your CRAB jobs just using the following (remember to set the correct environment for CRAB and also to create a proxy to the GRID using your GRID certificate!).
```
python multicrab.py
```


### Production checks while executing
As CRAB and Grafana are not very *helpful* sometimes, a specific Python tool, checkcrab.py, has been developed to cope with different (most probably not all) CRAB errors and problems and to provide a more credible status report than Grafana usually does.

Once you have your production ongoing, you can check its status using this tool as follows (inside the *test* folder).
```
python checkcrab.py WORKAREA NCORES
```
WORKAREA refers to the work area of your crab tasks (created as *analysis_* + *name*, where *name* is the homonymous variable used in multicrab) and is a compulsory argument, whereas NCORES refers to the same as before, and is optional.

NOTE: sometimes, inside each CRAB task, some jobs might fail due to not having enough time to completely execute. When running checkcrab, these jobs will be indeed detected and you will be asked whether you want to resubmit those jobs or not. However, at this moment, the script does not detect the error itself and it performs a "pure" resubmit. This will  be included in the future into the script, but to succesfully execute those tasks, you should do manually a resubmit with the --maxjobruntime option, setting for example --maxjobruntime 2000.

### Production checks after executing
Once your execution finishes, the script nanoAOD_checker.py allows you to check all produced datasets by first comparing the number of files that are produced with the number of files that are of MINIAOD, and second, by checking that the number of produced events is the same as the total number of events of the MINIAOD.

To use it, you should, in your storage site, establish any CMSSW-supported area, as well as the necessary CRAB & GRID environment and then copy the nanoAOD_checker.py script to the folder where all the production is stored (i.e. the folder named after the "name" variable in the multicrab.py script). Or, you can put it anywhere and then use the --folder (or -f) option to execute it. Then, simply do the following.
```
python nanoAOD_checker.py -n NCORES
```
This will start checking all the datasets in that folder. It is highly recommended to parallelise it, as it needs to loop over **all the files for each dataset**. The result is given at the end. You can set the debug option to true or false inside to get more or less logging.


