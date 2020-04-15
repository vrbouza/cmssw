import sys, os
import warnings as wr
import CentralSettings as cs
from multiprocessing import Pool

### Settings of the launch procedure: THEY ARE EXPECTED TO BE CHANGED DEPENDING ON THE EXECUTION

# This will be used to create the CRAB workarea, as well as the output folder
name = '2020_04_15_prodcheckout'

# The name will be append to the following path to set the output directory
#outputdir = '/store/user/rodrigvi/'
outputdir = '/store/user/rodrigvi/nanoAOD_production'

list_of_years = [
  #"2016",
  "2017",
  #"2018",
  ]

list_of_dataset_groups = [
    "All",   # All dataset groups
    #"DY",
    #"TT",
    #"WJets",
    #"WW",
    #"WZ",
    #"ZJ",
    #"SingleMuon",
    #"DoubleMuon",
    #"LowEGJet",
    #"HighEGJet",
]

if __name__ == '__main__':
    # Preventive errors.
    if   len(list_of_years) == 0:          raise RuntimeError("FATAL: no years chosen.")
    elif len(list_of_dataset_groups) == 0: raise RuntimeError("FATAL: no dataset group chosen.")

    # Getting task list; task := (dataset, dataset_group, year, crab_workarea = "prod_" + name, output_folder = outputdir + "/" + name)
    tasks = cs.GetAllTheDatasets(list_of_dataset_groups, list_of_years, "prod_" + name, outputdir + "/" + name)

    print "\n> We are going to launch the following {n} CRAB tasks in the workarea {w} that will drop their output in the folder {o}.".format(n = len(tasks), w = 'prod_' + name, o = outputdir + "/" + name + "/")
    for tsk in tasks: print "# Year {y}, dataset group {g}, dataset {d}.".format(d = tsk[0], g = tsk[1], y = tsk[2])
    print ""
    if not cs.confirm():
        print ""
        sys.exit()

    listoftaskswcreatedfolder = cs.CheckExistanceOfFolders(tasks)

    if len(listoftaskswcreatedfolder) != 0:
        print "\nWARNING: there are {n} task/s that we are going to run that were/was previously executed in the chosen workarea. The CRAB task submission will fail.".format(n = len(listoftaskswcreatedfolder))
        if not cs.confirm("Do you wish to send (just in case) a kill order to those tasks and afterwards erase the folders and continue with the execution? "):
            print ""
            sys.exit()
        wr.warn("WARNING: have in mind that some outputs MIGHT already been created at the output destination.")

        if len(sys.argv) > 1:
            ncores = int(sys.argv[1])
            print "\n> Parallelisation of old CRAB task killing and erasing with", ncores, "cores"
            pool = Pool(ncores)
            pool.map(cs.KillAndErase, listoftaskswcreatedfolder)
            pool.close()
            pool.join()
            del pool
        else:
            for tsk in listoftaskswcreatedfolder: cs.KillAndErase(tsk)

    #sys.exit()

    if not os.path.isdir("./" + tsk[3]): os.system("mkdir ./" + tsk[3])

    if len(sys.argv) > 1:
        print "PARALLELISATION OF CRAB LAUNCHING TEMPORALY DISABLED"
        #ncores = int(sys.argv[1])
        #print "\n> Parallelisation of CRAB task launching with", ncores, "cores"
        #pool = Pool(ncores)
        #pool.map(cs.LaunchCRABTask, tasks)
        #pool.close()
        #pool.join()
        #del pool
    #else:
        #for tsk in tasks: cs.LaunchCRABTask(tsk)
    for tsk in tasks: cs.LaunchCRABTask(tsk)
