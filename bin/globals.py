
"""Declares global symbols for use in FLASH setup"""

# Symbols exported to main namespace for those who want it
__all__ = ["GVars", "SetupError"]

##################################################################
########### Internal Implementation 
##################################################################

import string, sys, os.path, types

# Different levels of verbosity
PPDEBUG = 10   # Used for debugging Config file preprocessing
PPWARN  = 20   # Prints preprocessor warnings as well 
DEBUG   = 30   # Use for "current value of variable x is this
INFO    = 40   # Information about what the code is doing (not needed on normal runs)
WARN    = 50   # WARNINGs to the user
IMPINFO = 60   # Important information which user should usually be told
ERROR   = 70   # Critical Error which must be flagged to user
ALWAYS  = 80   # Information which must always be displayed (SUCESS / FAILURE)

# Different grid interpolations
GRID_INTERP_MONOTONIC = "GRID_INTERP_MONOTONIC"
GRID_INTERP_NATIVE    = "GRID_INTERP_NATIVE"

# Different geometries
GRID_GEOM_CARTESIAN   = "GRID_GEOM_CARTESIAN"
GRID_GEOM_CYLINDRICAL = "GRID_GEOM_CYLINDRICAL"
GRID_GEOM_SPHERICAL   = "GRID_GEOM_SPHERICAL"
GRID_GEOM_POLAR       = "GRID_GEOM_POLAR"
GRID_GEOM_UNDEF       = "GRID_GEOM_UNDEF"

# This string is used many times in the code
# reduces chance of typos
SimSetupDirname         = 'Simulation/SimulationMain'
SimRootDir              = 'Simulation'
SetupUnitsFilename      = 'setup_units'
UnitsFilename           = 'Units'
SetupParamsFilename     = "setup_params"
SetupVarsFilename       = "setup_vars"
SetupDefinesFilename    = "setup_defines"
SetupCallFilename       = "setup_call"
SetupDatafilesFilename  = "setup_datafiles"
FlashDefinesFilename    = "Simulation.h"
SimIntToStrFilename     = 'Simulation_mapIntToStr.F90'
SimStrToIntFilename     = 'Simulation_mapStrToInt.F90'
SimParticlesVarFilename = 'Simulation_mapParticlesVar.F90'
BuildStampGenFilename   = 'make_bstamp'
SetupLibrariesFilename  = 'setup_libraries'
RPInitParamsFilename    = 'rp_initParameters.F90'
RPDefaultParFilename    = 'default.par'
SetupFlashUnitsFilename = 'setup_getFlashUnits.F90'
RenormGroupFilename     = 'Simulation_getRenormGroup.F90'
VarnameTypeFilename     = 'Simulation_getVarnameType.F90'
Pm3RPFilename           = 'amr_runtime_parameters'
SuccessFilename         = ".success"
SimParticleTypeFilename = 'Particles_specifyMethods.F90'
EosMapFilename          = 'eos_variableMap.F90'

# Names of template files. Usually located in bin directory
MakefileTemplate        = 'Makefile.tpl'
BuildStampTemplate      = 'BuildStamp.tpl'
FlashDefinesTemplate    = 'Simulation.h.tpl'
SimIntToStrTemplate     = 'MapIntStr.F90.tpl'
SimStrToIntTemplate     = 'MapStrInt.F90.tpl'
SimParticlesVarTemplate = 'getParticlesVarMap.F90.tpl'
SetupFlashUnitsTemplate = 'FlashUnits.F90.tpl'
MakefileStubTemplate    = 'MakefileStub.tpl'
RenormGroupTemplate     = 'GetRenormGroup.F90.tpl'
VarnameTypeTemplate     = 'GetVarnameType.F90.tpl'
Pm3RPTemplate           = 'amr_runtime_parameters.tpl'
SimParticleTypeTemplate = 'Particles_specifyMethods.F90.tpl'
EosMapTemplate          = 'eos_variableMap.F90.tpl'

# List of prefixes for files created by setup. files with following prefixes 
# should not be touched by noClobber
noClobberExceptionList = [ "setup_", FlashDefinesFilename, SimIntToStrFilename, SimStrToIntFilename, 
                      UnitsFilename, BuildStampGenFilename, RPInitParamsFilename, RenormGroupFilename ]

SHORTCUT_CHAR = "+"

# for processing internal and external libraries
COMPILERS = ('FFLAGS', 'CFLAGS', 'LFLAGS', 'LIB')
DEFLTFLAG = 'OPT'

# unit prefixes which setup knows about
gridPrefix = "Grid/GridMain"
gridChoices = ["UG","UGSpatialFirst","AMR","Samrai","Chombo"]
simulationPrefix = "Simulation/SimulationMain"

# unit keywords for load unit specific modules and enforcing
# tomlfile design. New units should be registered here.
UNIT_KEYWORDS = ["Eos","Hydro","Inlet","Outlet","Simulation","Grid",
                 "Driver","IncompNS","HeatAD","Multiphase","Hydro",
                 "SolidMechanics","ImBound","Gravity","RadTrans",
                 "Spacetime","TimeAdvance","IO","Particles",
                 "PhysicalConstants","Logfile","Timers","Profiler",
                 "Burn","Deleptonize","Orchestration","Heater","Stir"]

######## Class for SetupError Exception
class SetupError(Exception):
    def __init__(self, message):
        super().__init__(message+"\nERROR")

####### Pretty printing class
class IndentedOutput:
    def __init__(self, numSpaces=4, file=None,debuglevel=WARN):
        self.WRAP = 80
        self.indent = 0
        self.numSpaces = numSpaces
        self.debuglevel = debuglevel
        if not file: self.file = sys.stdout
        else: self.file = file            

    def setDebugLevel(self,level):
        self.debuglevel = level

    def setWrapCol(self,wrap):
        self.WRAP = wrap

    # Return a list of strings so that each string has length
    # <= wrap and consist of full words only. If one word is too
    # long it may be an entry in the list
    def __pretty(self,text,wrap):
        ans = [] # soln so far
        more = text
        while more:
              if len(more) <= wrap: 
                 ans.append(more)
                 break
              # find last space before or equal to wrap length
              a = more.rfind(" ",0,wrap+1)
              if a < 0: # first word is long
                 b = more.find(" ")
                 if b < 0: 
                    ans.append(more)
                    more = ""
                 else:
                    ans.append(more[:b])
                    more = more[b:]
              else: # a >= 0
                 ans.append(more[:a])
                 more = more[a:]
              # if anything left remove leading spaces
              more = more.strip()
        return ans

    def put(self, text, level=WARN): 
        if level < self.debuglevel: return # info not important
        for line in map(str.strip, str.split(text, '\n')):
            nsp = self.indent*self.numSpaces
            ans = self.__pretty(line,self.WRAP-nsp)
            if not ans: # printing an empty string
               self.file.write('\n')
               self.file.flush()
               continue
            for a in ans:
                self.file.write("%s%s\n"% (' '*nsp,a))
                self.file.flush()

    def pop(self, numLevs=1):
        """Remove numLevs indentation levels"""
        self.indent = self.indent - numLevs
        if self.indent < 0:
            raise SetupError('Popped indentation one time too many!')

    def push(self, numLevs=1):
        self.indent += numLevs

######### Handle setup variables

class SetupVarsClass:
    vars = {} # actually contains the variables and values

    def __init__(self):
        self.set("fixedBlockSize", True)    # fbs is true by default

    def getdict(self): # get all the values as a dictionary
        return self.vars

    def addunit(self,unitname,val=False): # add a with variable for toplevel unit of specified unit
        isupper = lambda x: x == x.upper()
        parts = unitname.split(os.sep)
        while parts and parts[0] and not isupper(parts[0][0]): del parts[0]
        if parts:
           self.set("with%s"%parts[0],val)

    # set key to value, value may be a quoted string
    def set(self,key,value):
        if not value: value = ""
        if not isinstance(value,str): # direct store
           self.vars[key]=value
           return
        # we have a string which may require type change
        if len(value) >= 2 and value[0] == value[-1] and value[0] in ["'",'"']:
           value = str(value[1:-1])
        elif value.lower() in ["true","false","yes","no"]: # make it boolean
           value = value.lower() in ["true","yes"]
        else:
           try:
              a = int(value)
           except: pass
           else: value = a
        self.vars[key]=value

    def get(self,key):
        return self.vars.get(key,"")

    # file is a file object opened for writing
    def printvars(self,file,ignoreprefix="with"):
        vi = list(self.vars.items())
        vi.sort()
        for (k,v) in vi:
           if ignoreprefix and k.startswith(ignoreprefix): continue
           file.write('%s --> %s %s\n' % (k,v,type(v)))
        file.write('\n')
        

######### Global Variables

class GVarsClass:
    """Stores Global variables (visible to most of the code). Also includes parsed version of
    command line options"""
    out          = IndentedOutput()  # pretty printer
    wrapcol      = None
    setupVars    = SetupVarsClass() # handles setup variables
    indexReorder = False  # reorder indices in unk or not

    portable     = 0      # copy or link files
    verbose      = WARN    # level of verbose printing
    report       = 0

    dimension    = 2     # dimensionality of problem trying to solve
    defines      = []      # PreProcessor symbols to inform compiler, names and values if given
    definesNames = []      # PreProcessor symbols to inform compiler, names only
    maxblocks    = None  
    makedisplay  = 1   # dont display lines as they are executed in makefile
    makefileext  = ""  # link to Makefile.h+Thisextension
    nxb          = None  # domain size
    nyb          = None
    nzb          = None
    build_site   = None 
    build_tau    = None
    build_os     = None
    auto         = 0
    simulationName = None
    objectDir    = 'object'
    buildFlag    = 'OPT' #other options are 'DEBUG' and 'TEST' and 'OPT'
    npg          = 0
    withUnits    = []
    withUnofficial   = []
    withLibraries = {}
    datafiles    = [] # list of unix patterns for files to be copied over to object directory
    parfile      = "" # list of unix patterns for files to be copied over to object directory
    noClobber    = 0
    killUnits    = {} # units to be removed after all conditions have been met (USER BEWARE)
    withoutUnits = {} # units not needed (unless really needed)
    withoutLibraries = {} # for completeness store all withoutLibrary options (now dont do anyting with it)
    unitsFile    = "" # name of Units file to pick the units from
    gridInterpolation = None
    gridGeometry = GRID_GEOM_UNDEF
    curvilinear  = 0
    strictParams = 0
    particleMethods = {}
    particlesUnitIncluded = False
    nonexistent = "NONEXISTENT"
    eosStaticList = ['PRES','DENS','EINT','TEMP','GAMC','GAME','ENER','VELX','VELY','VELZ','SUMY','YE','ENTR','PRES1','PRES2','PRES3','EINT1','EINT2','EINT3','TEMP1','TEMP2','TEMP3','E1','E2','E3','SELE','SRAD']
    strEos = ''.join([x.lower()+'|'+x.upper()+'|' for x in eosStaticList])[:-1]
    macroOnly = False
    tomlfile = None
    tomlDict = {}
    withUnitMods = False

    def init(self,flashHomeDir):
        self.flashHomeDir = flashHomeDir
        self.sourceDir = os.path.join(self.flashHomeDir, 'source')
        self.simulationsDir = os.path.join(self.sourceDir, SimSetupDirname)
        self.libDir = os.path.join(self.flashHomeDir, 'lib')
        self.binDir = os.path.join(self.flashHomeDir, 'bin')
        self.tauDir = os.path.join(self.flashHomeDir, 'tools', 'tau')
        cwd = os.getcwd()
        os.chdir(self.sourceDir) # go to source Dir
        # remove initial "./" from each dirname
        self.topUnitNames = [x[2:] for x in self.getTopUnitNames() ]
        os.chdir(cwd)

    def getTopUnitNames(self, path='.',ignoreUnits=["flashUtilities"]):
      """Return the list of all top level unit names """
     
      # path does not refer to a directory
      if not os.path.isdir(path): return []

      basename = os.path.basename(os.path.abspath(path))
      ch = basename[0]
      if (ch >= 'A') and (ch <= 'Z'): # Are we a unit
         return [path]
      # we are only an organizational directory

      ans = []
      for dir in sorted(os.listdir(path)):
          if dir[0] == "." or (dir in ignoreUnits): continue
          subdir = os.path.join(path,dir)
          if not os.path.isdir(subdir): continue
          ans.extend(self.getTopUnitNames(subdir,ignoreUnits))
      return ans


############################################################################# 
##################### initialization code   ################################# 
############################################################################# 

GVars = GVarsClass()

