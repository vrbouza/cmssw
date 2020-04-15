import os, sys, argparse, json, copy
#import scipy.stats as stats
import subprocess as sp
import numpy as np
import ROOT as r
from multiprocessing import Pool, Manager

#### Directories with settings
all_data_dataset_groups = ["SingleMuon", "DoubleMuon", "SingleElectron", "DoubleEG", "MuonEG", "LowEGJet", "HighEGJet"]

r.gROOT.SetBatch(True);

debug = True


def get_datasets(path):
    finaldict = {}
    listofds  = []
    totalds   = 0
    for year in next(os.walk(path))[1]:
        finaldict[year] = {}
        for dg in next(os.walk(path + "/" + year))[1]:

            if len(next(os.walk(path + "/" + year + "/" + dg))[1]) == 0:
                finaldict[year][dg] = "NODATASETS"
            else:
                finaldict[year][dg] = {}
                for pref in next(os.walk(path + "/" + year + "/" + dg))[1]:
                    finaldict[year][dg][pref] = {}
                    for d in next(os.walk(path + "/" + year + "/" + dg + "/" + pref))[1]:
                        suff = d.replace(pref + "_", "")
                        finaldict[year][dg][pref][suff] = {}
                        totalds += 1
                        listofds.append(year + "/" + dg + "/" + pref + "/" + suff)
    return finaldict, totalds, listofds


def get_real_dataset_info(fulldataset):
    #print fulldataset
    #print "/".join(fulldataset.split("/")[1:])

    miniaodtag = "MINIAOD" + "SIM" * (fulldataset.split("/")[1] not in all_data_dataset_groups)

    #output_nfiles  = json.loads(sp.check_output('dasgoclient -query="file dataset=/' + "/".join(fulldataset.split("/")[2:]) + '/{tag}" -json'.format(tag = miniaodtag), shell = True))
    #output_nevents = json.loads(sp.check_output('dasgoclient -query="file,run,lumi,events dataset=/' + "/".join(fulldataset.split("/")[2:]) + '/{tag}" -json'.format(tag = miniaodtag), shell = True))

    output_nfiles  = json.loads(sp.check_output('dasgoclient -query="file dataset=/' + "/".join(fulldataset.split("/")[2:]) + '/{tag}" -json'.format(tag = miniaodtag), shell = True))
    output_nevents = json.loads(sp.check_output('dasgoclient -query="file,run,lumi,events dataset=/' + "/".join(fulldataset.split("/")[2:]) + '/{tag}" -json'.format(tag = miniaodtag), shell = True))

    nfiles = len(output_nfiles)


    totalevs = 0
    for f in output_nevents:
        for run in f[u'events']:
            totalevs += sum(run[u'number'])

    if debug:
        print "\n======================= DAS output"
        print "nfiles:",  nfiles
        print "nevents:", totalevs
        #print "output_nfiles:",  output_nfiles
        #print "output_nevents:", output_nevents
        print "=======================\n"

    #sys.exit()
    return nfiles, totalevs


def get_produced_dataset_info(fulldataset, folder, datadict, ncores = 1):
    yr = fulldataset.split("/")[0]; dg = fulldataset.split("/")[1]; pref = fulldataset.split("/")[2];
    suff = fulldataset.split("/")[3];
    dspath = "/".join([yr, dg, pref, pref + "_" + suff])
    subprods = os.listdir(folder + "/" + dspath)
    if debug: print "Initiating produced dataset checking"
    for subprod in subprods:
        datadict[yr][dg][pref][suff]["subprod_" + subprod] = {}

        if debug: print "\n====== Checking subprod:", subprod

        tmpnfiles  = 0
        tmpnevents = 0
        for subdir in next(os.walk(folder + "/" + dspath + "/" + subprod))[1]:
            if debug: print "### Checking subdir:", subdir
            tmppath = folder + "/" + dspath + "/" + subprod + "/" + subdir

            #tmplistoffiles = [get_tree_entries(tmppath + "/" + rootfile) for rootfile in os.listdir(tmppath) if (rootfile[-5:] == ".root")]

            tmptsks = [tmppath + "/" + rootfile for rootfile in os.listdir(tmppath) if (rootfile[-5:] == ".root")]

            pool = Pool(ncores)
            tmplistoffiles = pool.map(get_tree_entries, tmptsks)
            pool.close()
            pool.join()
            del pool

            tmpnevents += sum(tmplistoffiles)
            tmpnfiles  += len(tmplistoffiles)
        #print datadict
        datadict[yr][dg][pref][suff]["subprod_" + subprod]["nfiles"]  = int(tmpnfiles)
        datadict[yr][dg][pref][suff]["subprod_" + subprod]["nevents"] = int(tmpnevents)

        if debug: print "Detected nfiles:", tmpnfiles, ", detected nevents:", tmpnevents

    #print datadict[yr][dg][pref][suff]
    return


def get_tree_entries(rootfile):
    if debug: print "# Checking file:", rootfile
    tmpfile = r.TFile(rootfile, "READ")
    entries = tmpfile.Events.GetEntries()
    tmpfile.Close(); del tmpfile;
    return entries


def analyse_dataset(fulldataset, datadict):
    #### possible enhancements:
    # check per file rather than in total nevs
    # check file size and compare with the ones in GRID
    yr = fulldataset.split("/")[0]; dg = fulldataset.split("/")[1]; pref = fulldataset.split("/")[2];
    suff = fulldataset.split("/")[3];

    # Get subproductions
    subprods = []
    for key in datadict[yr][dg][pref][suff]:
        if key[:7] == "subprod": subprods.append(key)

    for subprod in subprods:
        equalnfiles  = False
        equalnevents = False
        # ...first, check if the number of files is the same
        if datadict[yr][dg][pref][suff][subprod]["nfiles"] == datadict[yr][dg][pref][suff]["das_nfiles"]:
            equalnfiles = True

        # ...and secondly, check the number of events
        if datadict[yr][dg][pref][suff][subprod]["nevents"] == datadict[yr][dg][pref][suff]["das_nevents"]:
            equalnevents = True

        if equalnfiles and equalnevents:
            datadict[yr][dg][pref][suff][subprod]["check_result"] = "CORRECT"
        else:
            datadict[yr][dg][pref][suff][subprod]["check_result"] = "WRONG:" + ("nfiles" if not equalnfiles else "nevents")

    if len(subprods) == 1:
        datadict[yr][dg][pref][suff]["check_result"] = datadict[yr][dg][pref][suff][subprods[0]]["check_result"]
    else:
        datadict[yr][dg][pref][suff]["check_result"] = "MULTIPLE|"
        tmpallok = True
        lastsubprod = get_latest_subprod(subprods)[0]
        if datadict[yr][dg][pref][suff][lastsubprod]["check_result"] != "CORRECT":
            datadict[yr][dg][pref][suff]["check_result"] += datadict[yr][dg][pref][suff][lastsubprod]["check_result"]
        else:
            datadict[yr][dg][pref][suff]["check_result"] += "CORRECT"
    return


def get_latest_subprod(listofsubprods):
    result = get_ordered_list(8, listofsubprods) # year
    #print "jeje"
    #print result
    if len(result) > 1:
        result = get_ordered_list(10, result) # month
        if len(result) > 1:
            result = get_ordered_list(12, result) # day
            if len(result) > 1:
                result = get_ordered_list(15, result) # hour
                if len(result) > 1:
                    result = get_ordered_list(17, result) # minute
                    if len(result) > 1:
                        result = get_ordered_list(19, result) # second
    return result


def get_ordered_list(index, thelist):
    #print thelist
    if len(thelist) <= 1: return thelist
    tmpdict = {}
    for el in thelist: tmpdict[el] = int(el[index:index + 2])
    listtoreturn =  sorted(thelist, key = lambda x: -tmpdict[x])
    if tmpdict[listtoreturn[0]] == tmpdict[listtoreturn[1]]:
        return [el for el in listtoreturn if tmpdict[el] == tmpdict[listtoreturn[0]]]
    else:
        return listtoreturn[:1]


def check_dataset(tsk, ncores = 1):
    fulldataset, datadict, folder = tsk

    if debug:
        print "\nInitiating check for:"
        print "Fulldataset:", fulldataset
        print "Datadict:", datadict
        print "Folder:", folder


    yr = fulldataset.split("/")[0]; dg = fulldataset.split("/")[1]; pref = fulldataset.split("/")[2];
    suff = fulldataset.split("/")[3];

    #### First, get DAS info
    das_nfiles, das_nevs = get_real_dataset_info(fulldataset)

    datadict[yr][dg][pref][suff]["das_nfiles"]  = das_nfiles
    datadict[yr][dg][pref][suff]["das_nevents"] = das_nevs


    #### Then, get your info
    get_produced_dataset_info(fulldataset, folder, datadict, ncores)

    #### Now, analyse
    analyse_dataset(fulldataset, datadict)
    print "> Checked /" + fulldataset

    #if "ZJToEEJ" in fulldataset: sys.exit()
    return


def print_check_results(thedict):
    print "\n> Check results:"
    wrongds = 0
    massivenods = 0
    for year in thedict:
        print "\n- YEAR: " + year
        for dg in thedict[year]:
            print "\n### DATASET GROUP: " + dg

            if datasetdict[year][dg] == "NODATASETS":
                print "There are no datasets in this dataset group folder. This might happen due to massive (i.e. in all jobs) execution errors."
                massivenods += 1
            else:
                for pref in thedict[year][dg]:
                    for suff in thedict[year][dg][pref]:
                        print thedict[year][dg][pref][suff]["check_result"] + " === /" + pref + "/" + suff
                        if "WRONG" in thedict[year][dg][pref][suff]["check_result"]: wrongds += 1

    if (wrongds + massivenods) == 0:
        print "> All datasets are correct!"
    else:
        if wrongds > 0:
            print "> A total of " + str(wrongds) + " datasets have not been correctly produced."
        if massivenods > 0:
            print "> A total of " + str(massivenods) + " dataset groups have not been correctly produced due to massive execution errors."

    if debug:
        print "Finaldict"
        print thedict

    return



if __name__=="__main__":
    parser = argparse.ArgumentParser(usage = "python nanoAOD_checker.py [options]", description = "Checker tool for the outputs of nanoAOD production (NOT postprocessing)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--folder',   '-f', metavar = 'folder', dest = "folder", help = 'folder where the nanoAOD datasets will be searched on. If not given, ./ is assumed.', required = False, default = "./")
    parser.add_argument('--ncores',   '-n', metavar = 'ncores', dest = "ncores", help = 'number of cores for parallelisation of checking.', required = False, default = 1)

    args   = parser.parse_args()
    folder = args.folder
    ncores = int(args.ncores)

    datasetdict, totalds, listofds = get_datasets(folder)

    #print datasetdict

    print "> The following " + str(totalds) + " datasets inside the corresponding dataset groups and inside the respective year/supergroup have been found inside the folder " + folder
    for year in datasetdict:
        print "\n- YEAR: " + year
        for dg in datasetdict[year]:
            print "\n### DATASET GROUP: " + dg
            if datasetdict[year][dg] == "NODATASETS":
                print "There are no datasets in this dataset group folder. This might happen due to massive (i.e. in all jobs) execution errors."
            else:
                for pref in datasetdict[year][dg]:
                    for d in datasetdict[year][dg][pref]:
                        print "/" + pref + "/" + d

    print "\n> Beginning check of produced nanoAOD"
    tasks = []

    #manageddict = Manager().dict(datasetdict)
    manageddict = datasetdict

    for el in listofds:
        tasks.append( (el, manageddict, folder) )

    #if ncores == 1:
        #print "- Sequential execution chosen"
        #for el in tasks: check_dataset(el)
    #else:
        #raise RuntimeError("ERROR: still not implemented")
        #print "- Parallelised execution chosen"
        #pool = Pool(ncores)
        #pool.map(check_dataset, tasks)
        #pool.close()
        #pool.join()
        #del pool

    for el in tasks: check_dataset(el, ncores)

    print_check_results(manageddict)
    print ""





