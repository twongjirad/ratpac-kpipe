import os,sys
import StringIO
import json

# tier 2 running utilities for KPipe rat-pac

def generate_condor_submtfile( jobid, outfilename, inputfiles, outputfiles, script, arguments, njobs=1 ):
    condorsub = StringIO.StringIO()
    print >> condorsub,"# Generated Tier 2 KPipe Condor submission file"
    print >> condorsub,"# --------------------------------------------"
    print >> condorsub,"Universe                = vanilla"
    print >> condorsub,"Environment             =\"HOSTNAME=$HOSTNAME\""
    print >> condorsub,"Requirements            = UidDomain == \"cmsaf.mit.edu\" && Arch == \"X86_64\" && HasFileTransfer"
    print >> condorsub,"Executable              = ",script
    print >> condorsub,"Arguments               = ",arguments
    print >> condorsub,"Input                   = /dev/null"
    print >> condorsub,"Output                  = condor_out/condor_job%d.out" % ( jobid )
    print >> condorsub,"Error                   = condor_err/condor_job%d.err" % ( jobid )
    print >> condorsub,"Log                     = condor.log"
    inputstr = ""
    iinput = 0
    for arg,inputfile in inputfiles.items():
        inputstr += inputfile
        if iinput!=len(inputfiles)-1:
            inputstr += ", "
        iinput += 1
    ioutput = 0
    outstr = ""
    for arg,outputfile in outputfiles.items():
        outstr += outputfile
        if ioutput!=len(outputfiles)-1:
            outstr += ", "
        ioutput += 1
    print >> condorsub,"transfer_input_files    = ",inputstr
    print >> condorsub,"transfer_output_files = ",outstr
    print >> condorsub,"should_transfer_files   = YES"
    print >> condorsub,"when_to_transfer_output = ON_EXIT"
    print >> condorsub,"Queue"
    print >> condorsub,"# --------------------------------------------"

    out = open( outfilename, 'w' )
    print >> out, condorsub.getvalue()

def generate_script_wrapper( outfilename, inputfiles, outputfiles, job_json ):
    script_wrapper = StringIO.StringIO()
    print script_wrapper >> "#!/bin/sh"
    print script_wrapper >> "Transfer input files"
    return

def parse_job_file( jsonfile ):
    f = open( jsonfile, 'r' )
    rawdata = f.read()
    parsed = json.loads( rawdata )
    return parsed

def make_arg_lists( job_json ):
    nargs = int(job_json["job"]["nargs"])
    start = int(job_json["job"]["startjob"])
    njobs = int(job_json["job"]["njobs"])

    inputlists = {}
    outputlists = {}
    arglists = {}

    for jobid in xrange(start,start+njobs):
        inputlists[jobid] = {}
        outputlists[jobid] = {}
        arglists[jobid] = ""
        for arg in job_json["job"]["args"]:
            if int(arg["njobids"])>0:
                jobarg = ( jobid*int(arg["njobids"]) )
                inputfile = arg[ "value" ] % jobarg
            else:
                inputfile = arg[ "value" ]
            if arg["type"]=="input":
                if arg["transfer"]=="True":
                    inputlists[jobid][ arg["argname"] ] = inputfile
                    arglists[jobid] += " "+os.path.basename(inputfile)
                else:
                    inputlists[jobid][ arg["argname"] ] = inputfile
                    arglists[jobid] += " "+inputfile
            else:
                outputlists[jobid][ arg["argname"] ] = inputfile
                arglists[jobid] += " "+inputfile
        for packagefile in job_json["job"]["packagefiles"]:
            inputlists[ jobid ][ os.path.basename(packagefile) ] = packagefile

    return inputlists,outputlists,arglists

def launch_jobs( job_json, inputlists, outputlists, arglists, check_for_out=True ):
    joblist = inputlists.keys()
    joblist.sort()
        
    script  = job_json["job"]["script"]
    os.system("mkdir -p condor_scripts")
    os.system("mkdir -p condor_out")
    os.system("mkdir -p condor_err")
    for jobid in joblist:
        condorfile = "condor_scripts/condor_submit_jobid%d.condor"%(jobid)
        script_wrapper = "condor_scripts/wrapper_job%d_%s"%(jobid, script)
        inputs = inputlists[ jobid ]
        outputs = outputlists[ jobid ]
        print inputs
        generate_condor_submtfile( jobid, condorfile, inputs, outputs, job_json["job"]["script"], arglists[jobid] )
        break

if __name__ == "__main__":
    example = "example_job.json"
    if len(sys.argv)==2:
        example = sys.argv[1]
    job_json = parse_job_file( example )
    inputlists, outputlists, arglists = make_arg_lists( job_json )
    launch_jobs( job_json, inputlists, outputlists, arglists )
    
