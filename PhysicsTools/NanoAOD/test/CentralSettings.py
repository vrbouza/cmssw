import os, sys, json, copy
### General settings
UnretrieveThreshold = 7
nanoAODv = 6

#### Directories with settings
all_data_dataset_groups = ["SingleMuon", "DoubleMuon", "SingleElectron", "DoubleEG", "MuonEG", "LowEGJet", "HighEGJet"]

#### Function definitions
def LaunchCRABTask(tsk):
    dataset, dataset_group, year, crab_workarea, output_folder = tsk
    from CRABAPI.RawCommand       import crabCommand
    from CRABClient.UserUtilities import config
    from CRABClient.JobType       import CMSSWConfig
    from multiprocessing          import Process

    def submit(cfg):
        res = crabCommand('submit', config = cfg)
        del res
        return

    print "\n# Launching CRAB task for year {y}, dataset group {g} and dataset {d}.".format(y = year, g = dataset_group, d = dataset)
    config = config()

    config.General.requestName = GetRequestName(dataset, year, crab_workarea)

    config.General.workArea        = crab_workarea
    config.General.transferLogs    = True
    config.General.transferOutputs = True

    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = "./nanoAODproduction_cfg.py"
    config.JobType.pyCfgParams = GetConfigFileParameters(dataset, dataset_group, year)

    #print config.JobType.pyCfgParams
    #sys.exit()

    config.JobType.maxMemoryMB = 2500
    config.JobType.allowUndistributedCMSSW = True

    config.Data.inputDataset = dataset
    config.Data.inputDBS     = 'global'
    config.Data.splitting    = 'FileBased'
    #config.Data.splitting    = 'Automatic'
    #config.Data.splitting    = "EventAwareLumiBased"
    config.Data.unitsPerJob  = 1
    #config.Data.unitsPerJob  = 2000
    #config.Data.totalUnits   = 6000

    #config.Data.publication  = False
    config.Data.publication  = True
    config.Data.allowNonValidInputDataset = True

    config.Data.outLFNDirBase    = output_folder + "/" + str(year) + "/" + dataset_group + "/"
    #config.Data.outputDatasetTag = "_".join(dataset.split("/")[1:-1])
    config.Data.outputDatasetTag = dataset.split("/")[2]

    #print dataset.split("/")[2]
    #sys.exit()


    config.Site.storageSite = 'T2_ES_IFCA'
    #config.Site.storageSite = 'T2_CH_CERN'
    #config.Site.blacklist   = ['T2_BR_SPRACE', 'T2_US_Wisconsin', 'T1_RU_JINR', 'T2_RU_JINR', 'T2_EE_Estonia']
    #config.Site.blacklist   = []
    #config.Site.ignoreGlobalBlacklist = True

    #sys.exit()

    p = Process(target = submit, args = (config,))
    p.start()
    p.join()

    #res = crabCommand('submit', config = config)
    #del config, res
    del config
    CMSSWConfig.configurationCache.clear() #### NOTE: this is done in order to allow the parallelised CRAB job submission. For further
                                           ## information, please check the code on [1], the commit of [2] and the discussion of [3].
                                           ## [1]: https://github.com/dmwm/CRABClient/blob/master/src/python/CRABClient/JobType/CMSSWConfig.py
                                           ## [2]: https://github.com/dmwm/CRABClient/commit/a50bfc2d1f32093b76ba80956ee6c5bd6d61259e
                                           ## [3]: https://github.com/dmwm/CRABClient/pull/4824
    return


def GetRequestName(d, y, wa):
    if len(wa + "_" + y + "_".join(d.split("/")[:-1])) < 100:
        return wa + "_" + y + "_".join(d.split("/")[1:-1])
    elif len((wa + "_" + y + "_" + d.split("/")[1] + "_XX" + d.split("/")[2][(len(d.split("/")[2]) - (99 - 4 - len(y) - len(wa) - len(d.split("/")[1]))):])) >= 100:
        return (wa + "_" + y + "_" + (d.split("/")[1])[:24] + "XX_XX" + d.split("/")[2][(len(d.split("/")[2]) - (99 - 6 - len(y) - len(wa) - 24)):])
    else:
        return (wa + "_" + y + "_" + d.split("/")[1] + "_XX" + d.split("/")[2][(len(d.split("/")[2]) - (99 - 4 - len(y) - len(wa) - len(d.split("/")[1]))):])



def confirm(message = "Do you wish to continue?"):
    """
    Ask user to enter y(es) or n(o) (case-insensitive).
    :return: True if the answer is Y.
    :rtype: bool
    """
    answer = ""
    while answer not in ["y", "n", "yes", "no"]:
        answer = raw_input(message + " [Y/N]\n").lower()
    return answer[0] == "y"


def CheckExistanceOfFolders(listoftasks):
    listoftaskswcreatedfolder = []
    for tsk in listoftasks:
        if os.path.isdir("./" + tsk[3] + "/crab_" + GetRequestName(tsk[0], tsk[2], tsk[3])):
            listoftaskswcreatedfolder.append(tsk)
    # task := (dataset, dataset_group, year, crab_workarea, output_folder = outputdir + "/" + name)
    return listoftaskswcreatedfolder


def KillAndErase(tsk):
    print "### Task with folder", "crab_" + GetRequestName(tsk[0], tsk[2], tsk[3])
    print "# Killing..."
    os.system("crab kill -d ./{wa}/{fl}".format(wa = tsk[3], fl = "crab_" + GetRequestName(tsk[0], tsk[2], tsk[3])))
    print "# Erasing..."
    os.system("rm -rf ./{wa}/{fl}".format(wa = tsk[3], fl = "crab_" + GetRequestName(tsk[0], tsk[2], tsk[3])))
    return


def GetAllTheDatasets(listofdgs, listofyears, crab_workarea, output_folder):
    # task := (dataset, dataset_group, year, crab_workarea, output_folder = outputdir + "/" + name)
    listofdatasets = []
    isdata = False

    if ["All"] == listofdgs:
        for year in listofyears:
            tmpdict = GetDictionaryFromJSON(year, True)
            for dg in tmpdict:
                for dataset in tmpdict[dg]:
                    if str(dataset[0]) == "#": continue
                    listofdatasets.append( (str(dataset), dg, year, crab_workarea, output_folder) )
            tmpdict = GetDictionaryFromJSON(year, False)
            for dg in tmpdict:
                for dataset in tmpdict[dg]:
                    if str(dataset[0]) == "#": continue
                    listofdatasets.append( (str(dataset), dg, year, crab_workarea, output_folder) )
    else:
        for dg in listofdgs:
            isdata = IsItData(dg)
            for year in listofyears:
                tmpdict = GetDictionaryFromJSON(year, isdata)
                if dg in tmpdict:
                    for dataset in tmpdict[dg]:
                        if str(dataset[0]) == "#": continue
                        listofdatasets.append( (str(dataset), dg, year, crab_workarea, output_folder) )
    return listofdatasets


def GetDictionaryFromJSON(year, isdata):
    if year.isdigit(): f = ("data" if isdata else "mc") + year + ".json"
    else:              f = year + ".json"
    with open("datasets_MINIAOD/" + f , "r") as jfile:
        return copy.deepcopy(json.load(jfile))


def IsItData(dg):
    for datadg in all_data_dataset_groups:
        if datadg in dg:
            return True
    return False


def GetConfigFileParameters(dataset, dg, y):
    pars = []
    cfgdict = {}
    if y.isdigit():
        y_ = ("data" if IsItData(dg) else "mc") + y
    else:
        y_ = y

    with open("datasets_MINIAOD/nanoAODv{v}_cfg.json".format(v = nanoAODv), "r") as jfile:
        cfgdict = copy.deepcopy(json.load(jfile))


    # Getting datinness
    if IsItData(dg):
        pars.append("IsData=True")
    else:
        pars.append("IsData=False")

    if y_ not in cfgdict:
        raise RuntimeError("ERROR: the information for the year/supergroup {y} could not be found in the JSON configuration file for nanoAOD version {v} while trying to get the parameters for dataset {d} inside dataset group {g}.".format(y = y_, v = nanoAODv, d = dataset, g = dg))

    # Getting era
    if isinstance(cfgdict[y_]["era"], dict):
        if "v1" in cfgdict[y_]["era"] and int(y) == 2017:
            if "RunIIFall17MiniAODv2" in dataset and not IsItData(dg):
                pars.append(str("Era=" + cfgdict[y_]["era"]["v2"]))
            else:
                pars.append(str("Era=" + cfgdict[y_]["era"]["v1"]))
    else:
        pars.append(str("Era=" + cfgdict[y_]["era"]))
    #print pars
    #sys.exit()
    # Getting global tag
    if not isinstance(cfgdict[y_]["global_tag"], basestring):
        raise RuntimeError("ERROR: attempted to obtain global tag for dataset {d} inside dataset group {g} for year/supergroup {y} with nanoAOD version {v}, but there are multiple ones implemented in the JSON config. file, and this is not supported in the code.".format(y = y_, v = nanoAODv, d = dataset, g = dg))
    else:
        pars.append(str("GlobalTag=" + cfgdict[y_]["global_tag"]))

    return pars
