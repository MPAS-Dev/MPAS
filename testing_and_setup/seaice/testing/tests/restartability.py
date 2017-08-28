#!/usr/bin/env python

import os, shutil
from compare_mpas_files import compare_files
from testing_utils import *

#-------------------------------------------------------------------------

def restartability(mpasDevelopmentDir, domainsDir, domain, configuration, check):

    testDir = "restartability.%s.%s" %(configuration,domain)

    # make a test directory
    create_test_directory(testDir)

    title = "Test: Restartability, Configuration: %s, Domain: %s" %(configuration,domain)

    print_colour(title, "title")

    logfile = open("log_test.txt","w")
    logfile.write(title)

    # base run
    nProcs = 16

    nmlChanges = {"cice_model": {"config_run_duration":'24:00:00'}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":"24:00:00"}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("base", mpasDevelopmentDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile) != 0):
        os.chdir("..")
        return 1

    # first restart run
    nProcs = 16

    nmlChanges = {"cice_model": {"config_run_duration":'12:00:00'}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":"12:00:00"}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("restart", mpasDevelopmentDir, domainsDir, domain, configuration, nmlChanges, streamChanges, nProcs, logfile) != 0):
        os.chdir("..")
        return 1

    # restart
    nProcs = 32

    nmlChanges = {"cice_model": {"config_start_time":"file"},
                  "restart":    {"config_do_restart":True}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}

    streamChanges = []

    if (restart_model("restart", nmlChanges, streamChanges, nProcs, logfile) != 0):
        os.chdir("..")
        return 1


    # compare
    restart_file = "restart.2000-01-02_00.00.00.nc"

    file1 = "./base/restarts/%s" %(restart_file)
    file2 = "./restart/restarts/%s" %(restart_file)

    ignoreVarname = ["cellsOnCell","verticesOnCell","edgesOnEdge","edgesOnCell"]
    if (check):
        ignoreVarname.append("testArrayParallelism")
        ignoreVarname.append("testArrayReproducibility")

    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile,ignoreVarname)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "restartability")

    logfile.close()

    os.chdir("..")

    return failed

#-------------------------------------------------------------------------
