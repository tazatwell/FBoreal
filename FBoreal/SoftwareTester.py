# DLW Feb 2018: inspired by PySP test software
import os
from os.path import join, dirname, abspath
import time
import subprocess
import difflib
import filecmp
import shutil
import unittest
import shlex
import pandas
import filecmp

rebasingAll = False  # just hack to rebase individual tests
progname = "SimApril2018.py"
csvname = "SimReport.csv" # basis for regression test
thisdir = dirname(abspath(__file__))
baselineDir = join(thisdir, "baselines")
instanceDir = join((dirname(thisdir)), "Instances")

class TestExamples(unittest.TestCase):

    def setUp(self):
        shutil.copy("../unix_FBPfunc5NODEBUG.so", "./FBPfunc5NODEBUG.so")
        self._tempfiles = []
        class_name, test_name = self.id().split('.')[-2:]
        self._csvbaselinename = join(baselineDir,
                                     class_name+"."+test_name+".csv")

    def _run_cmd(self, cmd):
        class_name, test_name = self.id().split('.')[-2:]
        print("%s.%s: Testing command: %s" % (class_name,
                                              test_name,
                                              str(cmd)))
        args = shlex.split(cmd)
        
        outname = os.path.join(thisdir,
                               class_name+"."+test_name+".out")
        self._tempfiles.append(outname)
        with open(outname, "w") as f:
            subprocess.check_call(args,
                                  stdout=f,
                                  stderr=subprocess.STDOUT)

    def _do_one(self, ifolder, fpref, adict):
        """
        Copy a forest file then run the command for the test.
        inputs:
            ifolder: string with input folder name
            fpref: string with prefix for Forest.asc
            adict: options dictionary
        """
        if fpref != "":
            shutil.copyfile(join(ifolder,fpref+"Forest.asc"), join(ifolder,"Forest.asc"))
        cmd = "python "+ progname
        for opt in adict:
            if adict[opt] is not None:
                cmd += " --"+opt + " " + adict[opt]
            else:
                cmd += " --"+opt
        self._run_cmd(cmd)

    def _rebase(self, fromdir):
        ff = join(fromdir, csvname)
        tt = self._csvbaselinename
        print ("rebase",ff, tt)
        shutil.copy(ff, tt)

    def _check_baseline(self, tmpdir):
        myguy = join(tmpdir, csvname)
        dfme = pandas.read_csv(myguy)
        
        baseguy = self._csvbaselinename
        dfbase = pandas.read_csv(baseguy)

        lastcolme = dfme.iloc[:,-1]
        lastcolbase = dfbase.iloc[:,-1]

        self.assertTrue(lastcolme.equals(lastcolbase))
               
    def _cleanup(self):
        for fname in self._tempfiles:
            if os.path.isfile(fname):
                try:
                    os.remove(fname)
                except OSError:
                    pass
            else:
                try:
                    shutil.rmtree(fname)
                except:
                    pass
                
        self._tempfiles = []

    def test_just_help(self):
        # just a smoke test for help
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"help": None}
        self._do_one("", "", adict)
        self._cleanup()
        
    def test_9CellsDet1k60min(self):
        # this one comes from 9CellsHom unlike most others
        ifolder = join(instanceDir, "9CellsHom")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "1",
                 "weather": "constant",
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0",
                 "Fire-Period-Length": "60"
        }
        self._do_one(ifolder, "", adict)
        if rebasingAll:
            self._rebase(tmpdir)
        self._check_baseline(tmpdir)
        self._cleanup()
        
    def test_9CellsRandConstWind1k60min(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "10", 
                 "weather": "constant",
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "60",
                 "seed": "1134"
        }
        self._do_one(ifolder, "1k", adict)
        if rebasingAll:
            self._rebase(tmpdir)
        self._check_baseline(tmpdir)
        self._cleanup()

    def test_9CellsFlipWind1k60min(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "10", 
                 "weather": "rows",  
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "60",
                 "seed": "1134"
        }
        self._do_one(ifolder, "1k", adict)
        if rebasingAll:
            self._rebase(tmpdir)
        self._check_baseline(tmpdir)
        self._cleanup()

    def test_9CellsSpotting(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "10", 
                 "weather": "rows",  
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "60",
                 "spotting-parameter-data-file-name": "SpotSample.json",
                 "seed": "1134"
        }
        self._do_one(ifolder, "1k", adict)
        if rebasingAll:
            self._rebase(tmpdir)
        self._check_baseline(tmpdir)
        self._cleanup()

    def test_9CellsScenarios(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "10", 
                 "weather": "rows",  
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "60",
                 "seed": "1134",
                 "scenarios": None,
                 "output-grid": None
        }
        self._do_one(ifolder, "1k", adict)
        # special test
        basefile = join(baselineDir, class_name+"."+test_name+".dat")
        testfile = join(join(tmpdir, "Scenarios"), "Scenario1.dat")
        
        if rebasingAll:
            shutil.copy(testfile, basefile)

        self.assertTrue(filecmp.cmp(testfile, basefile, shallow=False))
        self._cleanup()

    def test_9CellsPlots(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "1", 
                 "weather": "rows",  
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "60",
                 "seed": "1134",
                 "plot" : None
        }
        self._do_one(ifolder, "1k", adict)
        # special test
        basefile = join(baselineDir, class_name+"."+test_name+".png")
        testfile = join(join(join(tmpdir, "Plots"), "Plots1"), "forest0001.png")
        
        if rebasingAll:
            shutil.copy(testfile, basefile)

        self.assertTrue(filecmp.cmp(testfile, basefile, shallow=True))
        self._cleanup()

    #==== start one minute fire interval ====
    def test_9CellsFlipWind20by201min(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "10", 
                 "weather": "rows",  
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "1",
                 "seed": "1134"
        }
        self._do_one(ifolder, "20by20", adict)
        if rebasingAll:
            self._rebase(tmpdir)
        self._check_baseline(tmpdir)
        self._cleanup()

    def test_9CellsFlipWind10by101min(self):
        ifolder = join(instanceDir, "9CellsDLW")
        class_name, test_name = self.id().split('.')[-2:]
        tmpdir = os.path.join(thisdir, class_name+"_"+test_name)
        self._tempfiles.append(tmpdir)
        # use None and the argument for flags such as ignitions 
        adict = {"input-instance-folder": ifolder,\
                 "ignitions": None,
                 "nsims": "10", 
                 "weather": "rows",  
                 "sim-years": "1",
                 "ROS-Threshold": "1e-3",
                 "output-folder": tmpdir,
                 "ROS-CV": "0.13",
                 "Fire-Period-Length": "1",
                 "seed": "1134"
        }
        self._do_one(ifolder, "10by10", adict)
        if rebasingAll:
            self._rebase(tmpdir)
        self._check_baseline(tmpdir)
        self._cleanup()

if __name__ == '__main__':
    unittest.main()
