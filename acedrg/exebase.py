#!/usr/bin/env  ccp4-python
# Python script
#
#
#     Copyright (C) 2014 --- 2019 Fei Long,  G. Murshudov
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Library.
#
#====================================================================
## The date of last modification: 21/07/2016
#

import os,os.path,sys
import platform
import glob,shutil
import re,string
from optparse import OptionParser 
import time
import math
import select
import random

if os.name != 'nt':
    import fcntl
import signal

#################################################   

class CExeCode :
    """ A generic abstract base class that is to be inheritted by other classes that warp 
        a executable code. Basically this class defines two main methods that characterize the
        procedures a job by an executable code, i.e.  forming a command line object and then
        run the code according to the command line object"""

    def __init__( self ):

        # set pathes
 
        self.exitCode = 0

        # test 
        self._setCmdLineAndFile()

    def _setCmdLineAndFile(self ):         # will be overrided by individual derived classes
        self._cmdline = 'ls '
        self.batch_name = None
        self._log_name = None
            
    def subExecute(self, mode = 0, err_str ="" ):

        if self._log_name :
            if os.name == 'nt':
                self.logModeWin(mode, err_str)
            else: 
                self.logModeLin(mode, err_str)
        else :
            self.interactiveMode()

    def logModeLin(self, mode = 0, err_str=""):
        """ execute self._cmdline and output to a log file (self._log_name)  
        in a nonblocking way"""

        logfile = open(self._log_name, 'w')
        logfile.write("\n============ PROCESS INFORMATION =============\n")

        # spawn a sub-process for a job and connect to its input/output(and error)
        # streams using pipes

        try : 
            import subprocess
            subProcess = subprocess.Popen(self._cmdline, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
            subProcess.stdin.close()
            outfile_sub  = subProcess.stdout
            outfd_sub    = outfile_sub.fileno()
        except ImportError:
            print "module 'subprocess' has not be found, you use a version of PYTHON below 2.4"
            print "Import module 'popen2' instead, the calculations will not be affected"
            import popen2
            subProcess = popen2.Popen4(self._cmdline)
            subProcess.tochild.close()
            outfile_sub  = subProcess.fromchild
            outfd_sub    = outfile_sub.fileno()
    
        self._pid  = subProcess.pid
        pid_c      = self._pid  + 1
        p_name = os.path.basename(self._cmdline.strip().split()[0])
        t__log_name = os.path.basename(self._log_name)
        
       
        log_size = 0 
        log_max  = 500000000
        outfile_eof = 0
        if not mode:
            self._nonBlockingFile(outfd_sub)
            while not outfile_eof :
                ready_to_read, ready_to_write, in_error = \
                     select.select([outfd_sub],[],[],60.0)
                if outfd_sub in ready_to_read:
                    out_quantum = outfile_sub.read()
                    out_quantum_bt = len(out_quantum)
                    log_size = log_size + out_quantum_bt
                    if log_size > log_max:
                        logfile.write("\nTHE SIZE OF LOG FILE, %s, EXCEEDED THE LIMIT AND PROCESS STOPPED\n"%t__log_name)
                        logfile.write("The process running %s is killed\n"%p_name)
                        logfile.close()
                        os.kill(self._pid, signal.SIGKILL)
                        os.kill(pid_c, signal.SIGKILL)
                        print "The process running %s is killed " %p_name
                        time.sleep(1.0)
                        # be safe stop all process
                        sys.exit(1)
                    else :
                        if out_quantum == '':
                            outfile_eof = 1
                        logfile.write(out_quantum)
                        logfile.flush()
        else :
            self._nonBlockingFile(outfd_sub)
            while not outfile_eof :
                ready_to_read, ready_to_write, in_error = \
                     select.select([outfd_sub],[],[],10)
                if outfd_sub in ready_to_read:
                    out_quantum = outfile_sub.read()
                    if out_quantum == '':
                        outfile_eof = 1
                    logfile.write(out_quantum)
                    logfile.flush()
                    # allProcInfo.append(out_quantum)

        self.exitCode = subProcess.wait()
      
        # if mode:
        #    for item in allProcInfo:
        #        logfile.write(item)

        logfile.write("=========END OF PROCESS INFORMATION\n ==========\n")
        logfile.close()

        outfile_sub.close()

        if self.exitCode:
            print "#-----------------------------------------------------------------#"
            print "The process stoped.\nCheck the associated log file '%s'\nfor the error information!"%self._log_name
            print "#-----------------------------------------------------------------#"
            if not err_str:
                # if normal runtime errors, stop the program(using log. no err_str)  
                sys.exit(1)
                
        return True

    def logModeWin(self, mode = 0, err_str=""):
        """ execute self._cmdline and output to a log file (self._log_name)  
        under Windows system"""

        logfile = open(self._log_name, 'w')
        logfile.write("\n============ PROCESS INFORMATION =============\n")

        # spawn a sub-process for a job and connect to its input/output(and error)
        # streams using pipes

        try :
            import subprocess
            subProcess = subprocess.Popen(self._cmdline, shell=True, stdout=logfile)
            #subProcess.stdin.close()
            #outfile_sub  = subProcess.stdout
            #outfd_sub    = outfile_sub.fileno()
        except ImportError:
            print "module 'subprocess' has not be found, you use a version of PYTHON below 2.4"
            print "Import module 'popen2' instead, the calculations will not be affected"
            print "This should not be happening. Exiting Now."
            sys.exit()
            #import popen2
            #subProcess = popen2.Popen4(self._cmdline)
            #subProcess.tochild.close()
            #outfile_sub  = subProcess.fromchild
            #outfd_sub    = outfile_sub.fileno()

        self._pid  = subProcess.pid
        pid_c      = self._pid  + 1
        p_name = os.path.basename(self._cmdline.strip().split()[0])
        t__log_name = os.path.basename(self._log_name)

        self.exitCode = subProcess.wait()

        logfile.write("=========END OF PROCESS INFORMATION\n ==========\n")
        logfile.close()

        #outfile_sub.close()

        if self.exitCode:
            print "#-----------------------------------------------------------------#"
            print "The process stoped.\nCheck the associated log file '%s'\nfor the error information!"%self._log_name
            print "#-----------------------------------------------------------------#"
            if not err_str:
                # if normal runtime errors, stop the program(using log. no err_str)
                sys.exit(1)

        return True
    
    def _nonBlockingFile(self,fd):
        f_flag = fcntl.fcntl(fd, fcntl.F_GETFL,0)
        fcntl.fcntl(fd, fcntl.F_SETFL, f_flag | os.O_NONBLOCK)

    def interactiveMode(self):
        """ execute self._cmdline  just like what is done in a shell """
        
        os.system(self._cmdline)


# end


