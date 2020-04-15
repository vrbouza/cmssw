import sys, os
import CentralSettings as cs
import subprocess as sp
from multiprocessing import Pool
from datetime import datetime


def CheckCRABTaskStatus(crabdirpath):
    print "> Checking status for CRAB directory path", crabdirpath + "..."
    statusoutput = sp.check_output("crab status -d {d}".format(d = crabdirpath), shell = True)
    serverstatus         = ""
    schedstatus          = ""

    nSubJobs             = 0
    nIdleSubJobs         = 0
    nRunningSubJobs      = 0
    nTransferringSubJobs = 0
    nFinishedSubJobs     = 0
    nFailedSubJobs       = 0

    for line in statusoutput.splitlines():
        if "Status on the CRAB server:" in line:
            serverstatus = line.replace("Status on the CRAB server:", "").replace(" ", "").replace("\t", "")
        if "Status on the scheduler:" in line:
            schedstatus  = line.replace("Status on the scheduler:", "").replace(" ", "").replace("\t", "")
        if schedstatus == "" and "Cannot retrieve the status_cache file" in line: schedstatus = "STATUSUNREACHABLE"
        if "Jobs status" in line:
            nSubJobs = int(line.split("(")[-1][:-1].split("/")[-1])
        if ("idle"         in line and "Warning" not in line): nIdleSubJobs         = int(line.split("(")[-1][:-1].split("/")[0])
        if ("running"      in line and "Warning" not in line): nRunningSubJobs      = int(line.split("(")[-1][:-1].split("/")[0])
        if ("transferring" in line and "Warning" not in line): nTransferringSubJobs = int(line.split("(")[-1][:-1].split("/")[0])
        if ("finished"     in line and "Warning" not in line): nFinishedSubJobs     = int(line.split("(")[-1][:-1].split("/")[0])
        if ("failed"       in line and "Warning" not in line
            and "step" not in line and "(" in line):           nFailedSubJobs       = int(line.split("(")[-1][:-1].split("/")[0])


    #print "  - Task with server status:   ", serverstatus + "."
    #print "  - Task with scheduler status:", schedstatus + "."

    return (crabdirpath, serverstatus, schedstatus, "total:{tot}-idle:{idl}-running:{rn}-transferring:{trf}-finished:{fin}-failed:{fal}".format(tot = nSubJobs, idl = nIdleSubJobs, rn = nRunningSubJobs, trf = nTransferringSubJobs, fin = nFinishedSubJobs, fal = nFailedSubJobs))


def RelaunchCRABTask(crabdirpath):
    options  = crabdirpath.split("/")[-1].replace("crab_", "")

    sample   = options.split("_")[0] if "nu_" not in options else "_".join(options.split("_")[:2])
    ageing   = options.split("_")[1] if "nu_" not in options else options.split("_")[2]
    rpcuse   = options.split("_")[2] if "nu_" not in options else options.split("_")[3]
    scenario = ""
    outdir   = ""
    if "noage" not in ageing: scenario = "_".join(options.split("_")[3:]) if "nu_" not in options else "_".join(options.split("_")[4:])

    useyoungseg = False
    if "youngseg" in scenario:
        useyoungseg = True
        scenario = scenario.replace("youngseg_", "")


    print "\n# Task for the sample", sample, "with ageing", ageing, "with rpc use", rpcuse, "and ageing scenario", scenario, "and with young segments (not aged)" if useyoungseg else "and with aged segments."

    print "# Importing CRAB output directory..."
    logfile = open(crabdirpath + "/crab.log", "r")

    for line in logfile.readlines():
        if "config.Data.outLFNDirBase" in line:
            outdir = line.replace("config.Data.outLFNDirBase", "").replace(" ", "").replace("=", "").replace("'", "").replace('"', "").replace("\n", "")
            break
    if outdir == "": raise RuntimeError("FATAL: no CRAB output directory could be obtained from crab log.")
    logfile.close()

    print "# CRAB output directory is", outdir
    print "# Erasing CRAB workarea..."
    os.system("rm -rf " + crabdirpath)

    print "# Relaunching task..."
    cs.LaunchCRABTask( (sample, ageing + "_" + rpcuse + "_youngseg" * useyoungseg,
                        scenario, crabdirpath.split("/")[0], outdir) )
    return


def ResubmitCRABTask(crabdirpath):
    print "> Resubmitting CRAB task with directory path", crabdirpath + "..."
    os.system("crab resubmit -d {d}".format(d = crabdirpath))

    return


def GetSubJobsFromJob(inputstring):
    nidle = 0; nrun = 0; ntransf = 0; nfin = 0; nfail = 0; ntot = 0

    tmplist = inputstring.split("-")
    for el in tmplist:
        if   "idle"         in el: nidle   = int(el.split(":")[1])
        elif "running"      in el: nrun    = int(el.split(":")[1])
        elif "transferring" in el: ntransf = int(el.split(":")[1])
        elif "finished"     in el: nfin    = int(el.split(":")[1])
        elif "failed"       in el: nfail   = int(el.split(":")[1])
        elif "total"        in el: ntot    = int(el.split(":")[1])

    return (nidle, nrun, ntransf, nfin, nfail, ntot)


def GetYMDandHMSfromDate(date):
    dateymd = date.split(" ")[0].split("-"); datehms = date.split(" ")[1].split(":")
    year = dateymd[0]; month  = dateymd[1]; day    = dateymd[2];
    hour = datehms[0]; minute = datehms[1]; second = datehms[2];

    return [year, month, day, hour, minute, second]


def IsDate1SomeMinutesBeforeDate2(date1, date2, minutes = 7):
    times1 = GetYMDandHMSfromDate(date1)
    times2 = GetYMDandHMSfromDate(date2)
    IndeedItIs = True
    TimeScale  = 0

    while (IndeedItIs == True and TimeScale < 5): # We ignore seconds
        if   times1[TimeScale] > times2[TimeScale]:                     IndeedItIs = False
        elif times1[TimeScale] == times2[TimeScale] and TimeScale == 4: IndeedItIs = False
        elif times1[TimeScale] < times2[TimeScale]  and TimeScale != 4: break
        elif (times1[TimeScale] < times2[TimeScale]
             and abs(times1[TimeScale] - times2[TimeScale]) < minutes): IndeedItIs = False

        TimeScale += 1

    return IndeedItIs




if __name__ == '__main__':
    if len(sys.argv) == 1: raise RuntimeError("FATAL: no folder given to check.")

    basedir = sys.argv[1]
    if basedir[-1] == "/": basedir = basedir[:-1]

    print "\n> We are going to check the folder", basedir, "to see if there are CRAB tasks that failed on submission."
    print "\n----> REMEMBER TO SET CRAB PROPERLY!!! <----"

    listofsubdirs      = os.listdir(basedir)
    listofcrabdirpaths = [basedir + "/" + subdir for subdir in listofsubdirs if "crab_" == subdir[:5]]

    statuslist = []

    ncores = 0
    if len(sys.argv) >= 3:
        ncores = int(sys.argv[2])
        print "\n> Parallelising with", ncores, "cores."
        pool = Pool(ncores)
        statuslist = pool.map(CheckCRABTaskStatus, listofcrabdirpaths)
        pool.close()
        pool.join()
        del pool
    else:
        for crabdirpath in listofcrabdirpaths:
            statuslist.append(CheckCRABTaskStatus(crabdirpath))

    print "\n> Reviewing task' status..."
    totaltasks      = float(len(statuslist))
    relaunchinglist = []
    submittedlist   = []
    completedlist   = []
    onlytransflist  = []
    withfailedlist  = []
    notfinnotfaillist = []
    fullyidlelist   = []
    nSubJobsDict    = {}
    justlaunchedunretrievedlist = []
    justlaunchedlist= []

    nidletotal = 0.; nruntotal = 0.; ntransftotal = 0.; nfintotal = 0.; nfailtotal = 0.; ntotaljobs = 0
    currentdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    for job in statuslist:
        nidle, nrun, ntransf, nfin, nfail, ntot = GetSubJobsFromJob(job[3])
        nidletotal += nidle; nruntotal += nrun; nfintotal += nfin; nfailtotal += nfail; ntotaljobs += ntot
        nSubJobsDict[job[0]] = {
            "nidle"   : nidle,
            "nrun"    : nrun,
            "ntransf" : ntransf,
            "nfin"    : nfin,
            "nfail"   : nfail,
            "ntot"    : ntot
        }
        if   job[1] == "SUBMITFAILED":
            crabdir = job[0].split("/")[-1]
            print "# Task of CRAB directory", crabdir, "with SUBMITFAILED server status: adding to relaunching list."
            relaunchinglist.append(job[0])
        elif job[1] == "NEWoncommandSUBMIT":
            crabdir = job[0].split("/")[-1]
            justlaunchedlist.append(job[0])
        elif job[1] == "SUBMITTED":
            if job[2] == "STATUSUNREACHABLE":
                crabdir = job[0].split("/")[-1]

                # Checking if the tasks were just launched
                logfile        = open(job[0] + "/crab.log", "r")
                submitdate = "";

                for line in logfile.readlines():
                    if "Executing command: 'submit'" in line:
                        submitdate = line.split("DEBUG ")[-1].split(".")[0]
                        break
                logfile.close(); del logfile

                justlaunched = IsDate1SomeMinutesBeforeDate2(submitdate, currentdate, cs.UnretrieveThreshold)

                if justlaunched:
                    print "# Task of CRAB directory", crabdir, "with SUBMITTED server status but with unretrievable server status. It seems that it just has been launched (less than {thrsh} minutes).".format(thrsh = cs.UnretrieveThreshold)
                    justlaunchedunretrievedlist.append(job[0])
                else:
                    print "# Task of CRAB directory", crabdir, "with SUBMITTED server status but with unretrievable server status: adding to relaunching list."
                    relaunchinglist.append(job[0])
            else:
                submittedlist.append(job[0])
                if   job[2] != "COMPLETED" and job[2] != "FAILED" and (nrun != 0 or ntransf != 0):
                    notfinnotfaillist.append(job[0])
                elif job[2] != "COMPLETED" and job[2] != "FAILED":
                    fullyidlelist.append(job[0])


        if job[2] == "COMPLETED":
            completedlist.append(job[0])

        if nidle == 0 and nrun == 0 and ntransf != 0 and nfail == 0: onlytransflist.append(job[0])
        if nfail != 0: withfailedlist.append(job[0])

    nnotfinishedjobs = nidletotal + nruntotal + nfailtotal


    if len(justlaunchedunretrievedlist) != 0:
        print "\n### WARNING: some tasks have unretrievable scheduler status, but they seem to have been just launched (in the previous {thrsh} minutes). The former usually happens as a consequence of the latter. These tasks are:".format(thrsh = cs.UnretrieveThreshold)
        for tsk in justlaunchedunretrievedlist: print "#", tsk
        print ""
        if cs.confirm("Do you wish, nevertheless, to relaunch these tasks? (If you answer 'no', they will be still taken into account in the report)"):
            relaunchinglist += justlaunchedunretrievedlist


    print "\n========== GLOBAL CRAB REPORT =========="
    print "# Total tasks:                          ", int(totaltasks)
    print "# Submitted tasks:                      ", str(len(submittedlist)) + "/" + str(int(totaltasks)),   "(%3.1f "%(len(submittedlist)/totaltasks * 100) + "%)"
    print "# Just launched tasks:                  ", str(len(justlaunchedlist)) + "/" + str(int(totaltasks)),   "(%3.1f "%(len(justlaunchedlist)/totaltasks * 100) + "%)"
    print "# Tasks tagged for relaunch:            ", str(len(relaunchinglist)) + "/" + str(int(totaltasks)), "(%3.1f "%(len(relaunchinglist)/totaltasks * 100) + "%)"
    if len(justlaunchedunretrievedlist) != 0: print "# Submitted tasks with unretriveable sched. status:", str(len(justlaunchedunretrievedlist)) + "/" + str(int(totaltasks)), "(%3.1f "%(len(justlaunchedunretrievedlist)/totaltasks * 100) + "%)"
    print "# Tasks with transferring (only) jobs:  ", str(len(onlytransflist)) + "/" + str(int(totaltasks)),  "(%3.1f "%(len(onlytransflist)/totaltasks * 100) + "%)"
    print "# Tasks with failed jobs:               ", str(len(withfailedlist)) + "/" + str(int(totaltasks)),  "(%3.1f "%(len(withfailedlist)/totaltasks * 100) + "%)"
    print "# Tasks in normal execution:            ", str(len(notfinnotfaillist)) + "/" + str(int(totaltasks)),  "(%3.1f "%(len(notfinnotfaillist)/totaltasks * 100) + "%)"
    print "# Tasks submitted, but completely idle: ", str(len(fullyidlelist)) + "/" + str(int(totaltasks)),  "(%3.1f "%(len(fullyidlelist)/totaltasks * 100) + "%)"
    print "# Completed tasks:                      ", str(len(completedlist)) + "/" + str(int(totaltasks)),   "(%3.1f "%(len(completedlist)/totaltasks * 100) + "%)"
    print ""
    print "# Total jobs:        ", int(ntotaljobs)
    if int(ntotaljobs) != 0:
        print "# Idle jobs:         ", str(int(nidletotal)) + "/" + str(int(ntotaljobs)), "(%3.1f "%(nidletotal/ntotaljobs * 100) + "% of total" + ((", %3.1f "%(nidletotal/nnotfinishedjobs * 100) + "% of not completed nor transferring") if (nnotfinishedjobs != 0) else "") + ")"
        print "# Running jobs:      ", str(int(nruntotal)) + "/" + str(int(ntotaljobs)), "(%3.1f "%(nruntotal/ntotaljobs * 100) + "% of total" + ((", %3.1f "%(nruntotal/nnotfinishedjobs * 100) + "% of not completed nor transferring") if (nnotfinishedjobs != 0) else "") + ")"
        print "# Transferring jobs: ", str(int(ntransftotal)) + "/" + str(int(ntotaljobs)), "(%3.1f "%(ntransftotal/ntotaljobs * 100) + "%)"
        print "# Finished jobs:     ", str(int(nfintotal)) + "/" + str(int(ntotaljobs)), "(%3.1f "%(nfintotal/ntotaljobs * 100) + "%)"
        print "# Failed jobs:       ", str(int(nfailtotal)) + "/" + str(int(ntotaljobs)), "(%3.1f "%(nfailtotal/ntotaljobs * 100) + "%)"

    print "\n========== DETAILED CRAB REPORT =========="
    print "> Tasks tagged for relaunch (not even submitted or submitted w/o error but with unknown and unretrievable scheduler status):"
    if len(relaunchinglist) == 0: print "# There is not any CRAB task in that situation!"
    else:
        print "### There are", len(relaunchinglist), "CRAB tasks in such situation."
        for job in relaunchinglist: print "#", job

    print "\n> Tasks submitted, but completely idle:"
    if len(fullyidlelist) == 0: print "# There is not any CRAB task in such situation!"
    else:
        print "### There are", len(fullyidlelist), "CRAB tasks in such situation."
        for job in fullyidlelist: print "#", job

    print "\n> Tasks in normal execution:"
    if len(notfinnotfaillist) == 0: print "# There is not any CRAB task in normal execution!"
    else:
        print "### There are", len(notfinnotfaillist), "CRAB tasks in normal execution."
        for job in notfinnotfaillist: print "#", job

    print "\n> Tasks with only transferring jobs:"
    if len(onlytransflist) == 0: print "# There is not any CRAB task with only transferring jobs!"
    else:
        print "### There are", len(onlytransflist), "CRAB tasks with only transferring jobs."
        for job in onlytransflist:
            print "### Task", job
            print "# Transferring jobs:", str(nSubJobsDict[job]["ntransf"]) + "/" + str(nSubJobsDict[job]["ntot"]), "(%3.1f "%(nSubJobsDict[job]["ntransf"]/float(nSubJobsDict[job]["ntot"]) * 100) + "%)"

    print "\n> Tasks with failed jobs:"
    if len(withfailedlist) == 0: print "# There is not any CRAB task with failed jobs!"
    else:
        print "### There are", len(withfailedlist), "CRAB tasks with failed jobs."
        for job in withfailedlist:
            print "### Task", job
            print "# Idle jobs:        ", str(nSubJobsDict[job]["nidle"]) + "/" + str(nSubJobsDict[job]["ntot"]), "(%3.1f "%(nSubJobsDict[job]["nidle"]/float(nSubJobsDict[job]["ntot"]) * 100) + "%)"
            print "# Running jobs:     ", str(nSubJobsDict[job]["nrun"]) + "/" + str(nSubJobsDict[job]["ntot"]), "(%3.1f "%(nSubJobsDict[job]["nrun"]/float(nSubJobsDict[job]["ntot"]) * 100) + "%)"
            print "# Transferring jobs:", str(nSubJobsDict[job]["ntransf"]) + "/" + str(nSubJobsDict[job]["ntot"]), "(%3.1f "%(nSubJobsDict[job]["ntransf"]/float(nSubJobsDict[job]["ntot"]) * 100) + "%)"
            print "# Finished jobs:    ", str(nSubJobsDict[job]["nfin"]) + "/" + str(nSubJobsDict[job]["ntot"]), "(%3.1f "%(nSubJobsDict[job]["nfin"]/float(nSubJobsDict[job]["ntot"]) * 100) + "%)"
            print "# Failed jobs:      ", str(nSubJobsDict[job]["nfail"]) + "/" + str(nSubJobsDict[job]["ntot"]), "(%3.1f "%(nSubJobsDict[job]["nfail"]/float(nSubJobsDict[job]["ntot"]) * 100) + "%)\n"

    if len(relaunchinglist) == 0 and len(withfailedlist) == 0:
        print "\n> No tasks marked for relaunching nor resubmitting. Exitting."
        sys.exit()

    print "\nThere are a total of", len(relaunchinglist), "tasks marked for relaunching and", len(withfailedlist), "tasks marked for resubmitting.\n"
    if not cs.confirm():
        print ""
        sys.exit()

    if len(relaunchinglist) != 0:
        print "\n>  Initiating relaunch..."
        if ncores != 0:
            pool = Pool(ncores)
            pool.map(RelaunchCRABTask, relaunchinglist)
            pool.close()
            pool.join()
            del pool
        else:
            for job in relaunchinglist: RelaunchCRABTask(job)
        print "\n> All tasks relaunched."

    if len(withfailedlist) != 0:
        print "\n> Initiating resubmit..."

        if ncores != 0:
            pool = Pool(ncores)
            pool.map(ResubmitCRABTask, withfailedlist)
            pool.close()
            pool.join()
            del pool
        else:
            for job in withfailedlist: ResubmitCRABTask(job)
        print "> All tasks resubmitted."
