
# Declares a class UnitList which contains a list of unit names
# together with methods to manipulate them

__all__ = ["UnitList"]

###########################################
import globals
from globals import *
from utils import *
import lazyFile
from tplUtils import *
import unitCfg
from unitCfg import * # for FlashUnit and UnitUnion classes
from rpUtils import * # for RPInfo Class
from varUtils import * # for VarInfo class

import re, os.path, types, string
try:
    from functools import cmp_to_key
except ImportError: pass

###########################################

class UnitList:

    def __init__(self, ignorePP=False):
       self.units = {} # dictionary mapping name to FlashUnit object
       self.ignorePP = ignorePP

    def checkEmpty(self):
        if self.units: # programming error
           GVars.out.put("Programming Error: self.units should be empty here",globals.ERROR)
           raise SetupError

    def removeUnits(self,unitnames):
    # Remove these units also
    # returns the number of units actually removed
       if type(unitnames)==type("string"): unitnames = [ unitnames]
       count = 0
       for unitname in unitnames:
           if not unitname: continue
           if unitname in self.units:
              count += 1
              del self.units[unitname]
       return count

    def removeSubTrees(self,unitnames):
    # remove these units and all their children
      if not unitnames: return # speed up common case
      rmunitnames = []
      for unitname in unitnames:
          if unitname.endswith(os.sep):
             rmunitnames.append(unitname[:-1])
          else: rmunitnames.append(unitname)
      # rmunitnames has all unitnames not ending with "/"
      rmlist = [] # names of units to be removed
      for unitname in rmunitnames:
          usep = unitname+os.sep
          for uname in list(self.units.keys()):
              if uname == unitname or uname.startswith(usep):
                 rmlist.append(uname)
      for uname in rmlist:
          del self.units[uname]

    def addUnits(self,unitnames, comment=""): 
       # Add these units also
       # returns the number of units actually added
       if isinstance(unitnames,str):
          unitnames = [ unitnames]
       count = 0
       for unitname in unitnames:
           if not unitname: continue
           if unitname in self.units: continue
           if count == 0 and comment:
              GVars.out.put(comment,globals.INFO)
              GVars.out.push()
           count += 1
           if comment: GVars.out.put(unitname,globals.IMPINFO)
           if unitname in GVars.topUnitNames:
              restrict = unitCfg.TOP_UNIT
           else:
              restrict = unitCfg.NO_RESTRICT
           self.units[unitname] = FlashUnit(unitname, restrict=restrict, ignorePP=self.ignorePP) # add unit

       if count > 0 and comment: GVars.out.pop()
       return count

    def hasUnit(self,unitname):
        return unitname in self.units

    def checkTopUnits(self):
        for topUnit in GVars.topUnitNames:
          initFile = os.path.basename(topUnit) + "_init.F90"
          if not os.path.isfile(os.path.join(topUnit, initFile)):
             GVars.out.put("\nWARNING: Missing 'init' File",globals.ERROR)
             GVars.out.push()
             GVars.out.put("Each API-level unit directory must directly contain an initialization file")
             GVars.out.put("whose name matches the name of that directory plus \"_init.F90\"")
             GVars.out.put("The directory \"%s\" does not contain a file \"%s\"" % (topUnit, initFile))
             GVars.out.pop()
             #raise SetupError("")  # uncomment this to cause an abort
                                    # (and change above "WARNING" to "ERROR")

    def addTopUnits(self):
        self.addUnits(GVars.topUnitNames)

    def checkSuggest(self):
       # list of units which have been suggested but not included
       badunits = []
       unames = list(self.units.keys())
       for uname in unames:
          for suglist in self.units[uname]['SUGGEST']:
              badlist = 1
              for sug in suglist:
                  if sug in unames: badlist = None
              if badlist: badunits.append((uname,suglist)) # append the alternatives
       if badunits: # some suggestion not taken
          GVars.out.put("\nIGNORED SUGGESTIONS",globals.IMPINFO)
          for (uname,suglist) in badunits:
              GVars.out.put("%s suggested one of the following units as well: " % uname,globals.IMPINFO)
              GVars.out.push()
              for sug in suglist: GVars.out.put(sug,globals.IMPINFO)
              GVars.out.pop()
          GVars.out.put("",globals.IMPINFO)

    def checkRequirements(self):
       for unit in list(self.units.values()):
           setsOfAlternatives = unit["REQUIRES"]  # will be a list of lists
           # each member of this list-within-a-list will represent
           # a series of alternative Units as specified in a Config
           # file by the syntax "REQUIRES A OR B OR C".
           for setOfAlternatives in setsOfAlternatives:
               for unitName in setOfAlternatives:
                   if unitName in self.units:
                       # We have one in the list of required units
                       # so don't bother checking the others
                       break
               else:
                   raise SetupError('%s REQUIRES %s, not included'%(unit.name, " or ".join(setOfAlternatives)))


    def checkExclusivity(self):
      GVars.out.put("Checking if Exclusive units are included",globals.DEBUG)
      #units contains everybody's parents
      for unit in list(self.units.values()):
        for group in unit['EXCLUSIVE']:
            #No two elements of group must be in units
            a=None
            for b in group:
                if b not in self.units: continue
                if a: raise SetupError('%s and %s are exclusive'%(a, b))
                a=b

    def checkConflicts(self):
      GVars.out.put("Checking for Conflicting units",globals.DEBUG)
      for unit in list(self.units.values()):
        for b in unit['CONFLICTS']:
            if b in self.units:
               raise SetupError("setup Error: requested unit %s CONFLICTS with %s" % (b,unit))
            
    def getUnitNamesFromFile(self,file):
      if not os.path.isfile(file):
         GVars.out.put('cannot access %s file'%file,globals.ERROR)
         GVars.out.push()
         GVars.out.put('* Either use the -auto option or',globals.ERROR)
         GVars.out.put('* specify the Units file using the -unitsfile=<filename> option',globals.ERROR)
         GVars.out.pop()
         raise SetupError('No Units file found', "NOUNITS")

      GVars.out.put('scanning %s file for included units'%file,globals.IMPINFO)
      ans = []
      for line in open(file).readlines():
        rawline = line
        if str.count(line, '#'):
            line=line[:str.find(line, '#')]
        line=str.strip(line)
        if not line: continue
        try:
            a,b=str.split(line)
            if a!='INCLUDE': raise SetupError
        except (ValueError, SetupError):
            raise SetupError('Bad syntax:\n%s'%rawline)
        ans.append(b)

      return ans
 
    def addUnitsFromFile(self):
        fname = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.UnitsFilename)
        self.addUnits(self.getUnitNamesFromFile(fname),"Scanning Units file")

    def getDefaultUnits(self): 
      # return names of default units of units already present
      return [unit["DEFAULT"] for unit in list(self.units.values())]

    def addDefaultUnits(self):
      # add default units of all units
      while self.addUnits(self.getDefaultUnits(),"") > 0: pass

    def addKernelUnits(self):
      #check to see which units are Kernel units
      #if we have a Kernel unit then we want to
      #include all units and directories under the
      #unit specified to be a Kernel
  
      ans = []
      GVars.out.push()
      GVars.out.put('*** KERNEL *** addKernelUnits *** START ***',globals.DEBUG)
      for unit in list(self.units.values()):
        if unit['KERNEL']: 
            GVars.out.put('Processing KERNEL for %s'%unit,globals.DEBUG)
            GVars.out.push()
            if not os.path.isdir(unit['KERNEL']):
                raise SetupError("Invalid KERNEL %s, not a directory, requested by %s.\n"%(unit['KERNEL'],unit))
            GVars.out.pop()
            ans.extend(self.recursiveGetDir(unit['KERNEL']))
      self.addUnits(ans)
      GVars.out.put('*** KERNEL *** addKernelUnits *** END ***',globals.DEBUG)
      GVars.out.pop()

    def addGridInterpolation(self):
        # Different interpolation units need to be added depending on which
        # Grid has been chosen. This is a somewhat hacky solution to that
        # problem. As far as I know this is the only instance of a command-
        # line option whose effect depends on which other Units have already
        # been included. -nttaylor

        # Note that no interpolation Unit needs to be added for UG
        ans = []
        if GVars.gridInterpolation == globals.GRID_INTERP_MONOTONIC:
            for unitName in list(self.units.keys()):
                if unitName.startswith("Grid/GridMain/AMR/Paramesh4"):
                    # Check if added this already
                    for unitName2 in list(self.units.keys()):
                        if unitName2.startswith("Grid/GridMain/AMR/Paramesh4/interpolation"):
                            break
                    else:
                        # add the required interpolation unit
                        ans.append("Grid/GridMain/AMR/Paramesh4/interpolation")
                    break
                elif unitName.startswith("Grid/GridMain/AMR/Paramesh2"):
                    # Check if added this already
                    for unitName2 in list(self.units.keys()):
                        if unitName2.startswith("Grid/GridMain/AMR/Paramesh2/monotonic"):
                            break
                    else:
                        # add the required interpolation unit
                        ans.append("Grid/GridMain/AMR/Paramesh2/monotonic")
                    break
        elif GVars.gridInterpolation == globals.GRID_INTERP_NATIVE:
            # This is the only other possiblity, and all it means is that we
            # don't have to add any special units, so there's nothing to do.
            pass

        self.addUnits(ans)

    def recursiveGetDir(self,path):
      GVars.out.put('Will walk %s ...'%path,globals.DEBUG)
      
      ans = []
      for root, dirs, files in os.walk(path):
          # Unfortunately these are not skipped automatically
          # but we can modify dirs in place for top-down traversal
          # of the directory structure to make sure that `os.walk`
          # properly excludes everything in `.*` directories
          for dir in dirs:
              if(dir[0] == '.'):
                  dirs.remove(dir)
          ans.append(root)
          
      # now ans contains the list of all directories inside path 
      # except "." directories and their children
      GVars.out.put('...and found %s by walking.'%ans,globals.DEBUG)
      return ans

    def getParentUnits(self):
        return [unit.getParent() for unit in list(self.units.values())]

    def addParentUnits(self):
        while self.addUnits(self.getParentUnits()) > 0: pass

    def removeParentUnits(self):
        self.removeUnits(self.getParentUnits())

    def getRequiredUnits(self): # list of unitnames to be added
      ans = []
      for unit in list(self.units.values()):
          setsOfAlternatives = unit["REQUIRES"]  # will be a list of lists
          # each member of this list-within-a-list will represent
          # a series of alternative Units as specified in a Config
          # file by the syntax "REQUIRES A OR B OR C". Only one of
          # these units need be appended to 'ans'
          for setOfAlternatives in setsOfAlternatives:
              for unitName in setOfAlternatives:
                  if unitName in self.units:
                      # we already have one in the list of required units
                      # so don't bother checking the others
                      break
              else:
                  # only if we didn't have any of the listed alternatives
                  # do we choose the first one by appending it to 'ans'
                  ans.append(setOfAlternatives[0])
      return ans


    def addRequiredUnits(self): # Satisfy requirements
      while self.addUnits(self.getRequiredUnits()) > 0: pass
                
    def writeUnitsfile(self):
      outfile = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.UnitsFilename)
      if os.path.isfile(outfile):
        GVars.out.push()    
        GVars.out.put('Backing up %s' % outfile,globals.INFO)
        os.rename(outfile, outfile+'.bak')
        GVars.out.pop()

      outfd = open(outfile, 'w')
      outfd.write('#Units file for %s generated by setup '\
                  '\n\n'%GVars.simulationName)

      # Dont print parent units (makes file concise)
      pUnits = self.getParentUnits()
      cUnits = [x for x in list(self.units.keys()) if x not in pUnits]
      cUnits.sort()
      for unitname in cUnits:
          outfd.write('INCLUDE %s\n'%unitname)
      outfd.close()
    

######################################### User calls methods below 

    def getList(self):
        tmpList = list(self.units.keys())
        tmpList.sort()
        return tmpList

    def getLinkOrder(self): 

      # Return true if child "unita" should be linked earlier
      def earlier(corder,childa,childb):
        # code comes here only if we have a non-trivial corder
        if childa in corder:
           try:
              return cmp(corder.index(childa),corder.index(childb))
           except: # unitb is not in the list
              return -1
        else: # unita not present in list
           if childb in corder:
              return 1
           else: # both not in corder
              return cmp(childa,childb)

      # default comparison between different children
      def defcmp(childa,childb):
          if childa == globals.SimRootDir: return 1
          if childb == globals.SimRootDir: return -1
          return cmp(childa,childb)

      # return 0 if equal, -1 if patha < pathb, 1 if patha > pathb
      def compare(apath,bpath):
          if apath == bpath: return 0
          # find youngest common ancestor
          common = []
          patha = apath[:] # copy is essential here
          pathb = bpath[:]
          while patha and pathb and (patha[0] == pathb[0]):
                common.append(patha[0])
                del patha[0]
                del pathb[0]
          parent = "/".join(common) # found youngest common parent
          # patha is a prefix of pathb
          if not patha: return -1
          # pathb is a prefix of patha
          if not pathb: return 1
          # Now need to compare childa and childb according to order specified by parent
          childa = patha[0]
          childb = pathb[0]
          # if trouble accessing parent's order (pick default ordering)
          try:
            corder = self.units[parent]['CHILDORDER']
            if corder: 
               return earlier(corder,childa,childb)
            else: return defcmp(childa,childb)
          except:
            return defcmp(childa,childb)
          
      # returns names of units present in a fixed order
      # split unitnames into list describing path
      pathlist = [unitname.split("/") for unitname in list(self.units.keys())]
      # First do a lexicographic sort as that should get most of elements in place
      # since this is optimized in python3, we have done most of the work without using
      # our fancy ordering
      pathlist.sort() 
      # now use our fancy sorting routine, so hopefully there are few calls to our 
      # comparator algorithm
      try:
          pathlist.sort( key = cmp_to_key(compare) )
      except NameError: # in case we did not import cmp_to_key successfully...
          pathlist.sort( cmp = compare )

      # combine paths back to strings
      tmpList = [ "/".join(x) for x in pathlist ] # contract paths to unitnames
      return tmpList

    def getConfigInfo(self,**kw): # Generate UnitUnion class based on given units
        cInfo = UnitUnion(list(self.units.values()),**kw)
        cInfo.unitNames = self.getList()
        return cInfo

    def getAllVariants(self,unitName):
      varList = []
      for var in self.units[unitName]['VARIANTS']:
        varUnit = unitName + "/" + var
        varList.append(varUnit)
      return varList

    def getRequestedVariants(self,unitName):
      varList = []
      for var in self.units[unitName]['VARIANTS']:
        varList.append(var)
      isNull = [ var.lower() == 'null' for var in self.units[unitName]['VARIANTS'] ]
      if (len(varList)==0 and any(isNull) ): # Huh? This can never be
          # true! What is the intention here? - KW
        varList.append('Null')
      return varList

    def recursiveGetDefs(self,sourceDir,targetUnit,requiringDir=None):
        """
        requiringDir: The directory containing the Config file that has the
        REQUIRES that got this function invoked. Only passed in when
        it is different from targetUnit. I.e., when targetUnit really
        comes from the target of an effective REQUIRES directive.
        """
        defsList = []
        # Since this list has been sorted, the children should be
        # checked in the proper order
        for unit in self.getLinkOrder():
            if(targetUnit == unit or
               targetUnit.startswith(unit+"/") or # targetUnit is unit+"/"+<something>
               unit.startswith(targetUnit+"/")    # unit is targetUnit+"/"+<something>
               ):
                isVariant = False
                if requiringDir is not None:
                    for varUnit in self.getAllVariants(requiringDir):
                        if unit == varUnit:
                            isVariant = True
                        elif unit.startswith(varUnit+"/"):
                            isVariant = True
                if not isVariant:
                    for varUnit in self.getAllVariants(targetUnit):
                        if unit == varUnit:
                            isVariant = True
                        elif unit.startswith(varUnit+"/"):
                            isVariant = True

                if(not isVariant):
                    unitDir = os.path.join(sourceDir,unit)
                    defsList += getDefs(unitDir)
                else:
                    GVars.out.push()
                    GVars.out.put("McDef files in %s skipped because of isVariant" % unit, globals.DEBUG)
                    GVars.out.pop()

        return defsList

    # A unit needs definitions from itself as well as any unit it REQUIRES.
    #
    # Directories that are ancestors of a target (REQUIRED) unit,
    # except those that are recognized as a variant or subvariant of the requiring unit,
    # are searched for McDef files.
    #
    # A REQUIRED unit itself,
    # except when it is recognized as a variant or subvariant of the requiring unit,
    # is searched for McDef files.
    #
    # Directories that are children or higher-degree descendants of a target (REQUIRED) unit,
    # except those that are recognized as a variant or subvariant of the current or the target unit,
    # are searched for McDef files.
    #
    # Ancestor directories of the current unit are searched for McDef files.
    # The current unit itself is searched for McDef files.
    # Directories that are children or higher-degree descendants of the current unit,
    # except those that are recognized as a variant or subvariant of the current unit,
    # are searched for McDef files.
    #
    # Subvariant of x: a directory that is either a variant of x or in
    # a directory subtree under x.
    # Cases that cannot happen are omitted from mention above. For example, it is taken
    # for granted that a variant of x cannot possibly be x itself or an ancestor of x.
    # The terms 'unit' and 'directory' are used interchangably, but it is assumed that
    # a 'unit' only matters if it is included in the 'unitList', i.e., the list of
    # directories that are included in the simulation configuration by the usual means.
    #
    # Note that the directories that are collected here and listed above include ONLY
    # those that exist in the unitList self, i.e., directories that are included in the
    # simulation configuration by the usual means of `-with-unit` or `-unit` setup
    # options, REQUIRES, REQUESTS, and DEFAULT directives (as far as they become
    # effective) in Config files, and equivalents thereof. These directories do NOT
    # include ones that are ONLY referenced by VARIANT directives in Config files.
    #
    # Also note that what is returned is not actually a list of directories, but a list
    # of two lists, each of them containing filenames of McDef files - for the most
    # part, these are the McDef files that are found by searching as described above.
    # Each McDef file is expected to occur at most once in a list.
    #
    # There is additional special handling for
    # - McDef files in the binDir: these should all appear at the beginning of the first list.
    # - McDef files in the simulation directory: if there are any, these should all
    #   appear as the contents of the second list.
    #
    # McDef files in directories that are not included in the simulation configuration
    # by the usual means, but are referenced by a VARIANT directive in an included
    # Config file, will be collected elsewhere! See generateVariants of the
    # macroProcessorHelper. The list of those "pure VARIANT" McDef files is to be
    # combined with the lists returned by the present function in such a way that, for
    # actual macro expansion, the expected rules of priority for handling same-name
    # macros from different McDef sources are obtained.
    def collectDefs(self,sourceDir,unitName,binDir,simDir):
        # We build, and shall return, two lists:
        defsList = [[],[]]
        # Variant-specific mc definitions (not collected here) should be
        # handled by the caller as if their priority was in between
        # the priorities of the first and the second returned list.

        GVars.out.push()
        GVars.out.put("Collect McDef files for unit %s in %s" % (unitName, sourceDir), globals.DEBUG)

        # get common defs from bin dir
        defsList[0] += getDefs(binDir)

        # get all defs from required units that are actually in the unitList
        requiredList = self.units[unitName]["REQUIRES"]
        # each member of this list-within-a-list will represent
        # a series of alternative Units as specified in a Config
        # file by the syntax "REQUIRES A OR B OR C".
        for setOfAlternatives in requiredList:
            # The following loop nest only descends into the **first** of
            # these required units that is actually present in self.units:
            for requiredName in setOfAlternatives:
                if requiredName in self.units:
                    for inifile in self.recursiveGetDefs(sourceDir,requiredName,unitName):
                        # prevent unintended overriding...
                        if inifile not in defsList[0]:
                            GVars.out.put("McDef file %s because unit %s REQUIRES %s" % (inifile, unitName, requiredName), globals.INFO)
                            defsList[0].append(inifile)
                    break

        # get defs for the unit itself
        for inifile in self.recursiveGetDefs(sourceDir,unitName):
            # prevent unintended overriding...
            if inifile not in defsList[0]:
                GVars.out.put("McDef file %s for unit %s itself" % (inifile, unitName), globals.INFO)
                defsList[0].append(inifile)
        # defs in simulation directory should overwrite others
        defsList[1] += getDefs(simDir)

        GVars.out.pop()
        return defsList

    def generateUnitsfile(self):
      # called only when '-auto' flag is passed
      GVars.out.put('generating default Units file',globals.IMPINFO)
      self.checkEmpty()

      # add the appropriate unit from "SimulationMain"
      self.addUnits(os.path.join(globals.SimSetupDirname, GVars.simulationName))
      # add units specified on command-line with "--with-unit"
      self.addUnits(GVars.withUnits)

      while 1:
        oldnames=self.getList()  # names of all units in this UnitList instance

        self.addParentUnits()
        self.addRequiredUnits()
        self.removeParentUnits()
        # start over if we added any required units
        if oldnames!=self.getList(): continue

        self.addDefaultUnits()
        self.addKernelUnits()

        self.addGridInterpolation()
        
        # Nothing new added so we are done
        if oldnames == self.getList(): 
           break
            
      # now kill units
      if GVars.killUnits:
         GVars.out.put('Killing subtrees: %s' % ",".join(list(GVars.killUnits.keys())),globals.INFO) # USER BEWARE
         self.removeSubTrees(list(GVars.killUnits.keys()))

      # Tell User about the list of units
      GVars.out.push()
      for unitname in self.getList():
          GVars.out.put(unitname,globals.IMPINFO)
      GVars.out.pop()

      self.writeUnitsfile()
      # If user had not used "-auto", we would not even be in this branch
      # of code, so now that the Units file is written, we clear out the
      # units dictionary again. Now we are in the same state regardless of
      # whether the "-auto" flag was passed or not. The "Units" file must
      # be read in by "addUnitsFromFile" in method "populate" below.
      self.units = {}

    def populate(self):
       # main stuff
       self.checkEmpty()
       simUnit = os.path.join(globals.SimSetupDirname, GVars.simulationName)

       self.addUnitsFromFile() # read Units file and add units found there
       if not GVars.auto: 
          # check if a Simulation is already present
          names = [x for x in self.getList() if x.startswith(globals.SimSetupDirname) ]
          try: # remove unit if present
            names.remove(simUnit)
          except ValueError:
            pass
          if names: 
             GVars.out.put("\nERROR: Multiple Simulation Units Found",globals.ERROR)
             GVars.out.push()
             GVars.out.put("You tried to setup %s as well as " % simUnit,globals.ERROR)
             GVars.out.put("use your previous (or specified) Units file.",globals.ERROR)
             GVars.out.put("Your units file already has %s. Your choices are \n"%names[0],globals.ERROR)
             GVars.out.put("* Use the -auto option and I will use my crystal ball",globals.ERROR)
             GVars.out.put("* Explicitly specify the units file using -unitsfile option",globals.ERROR)
             GVars.out.put("* Change the problem you want to setup",globals.ERROR)
             GVars.out.pop()
             raise SetupError("")
          else:
             self.addUnits(simUnit)

       GVars.out.put('checking for default sub-units in included units',globals.INFO)

       self.addDefaultUnits()
       self.addParentUnits()
       self.addKernelUnits()
    
       GVars.out.put('checking for exclusivity and conflicts',globals.INFO)
       self.checkExclusivity()
       self.checkConflicts()

       GVars.out.put('looking for paths of non-included units',globals.INFO)
       self.checkTopUnits()
       self.addTopUnits()

       self.addUnits(simUnit)

       GVars.out.put('checking requirements',globals.INFO)
       self.checkRequirements()

       self.checkSuggest()
       self.removeSubTrees(list(GVars.killUnits.keys())) # remove specified units and children (USER BEWARE no CHECK performed on these)

    def getRPInfo(self,max_plot_vars,**kw):
        rpInfo = RPInfo(**kw)
        for (unitname,unit) in list(self.units.items()):
            for (rpname,(rptype,rpvalue,rpconst,rprange)) in list(unit['PARAMETER'].items()):
                # try rpname, then rpname_parameter
                rpcomment = unit['D'].get(rpname,unit['D'].get("%s_parameter"%rpname,""))
                rpInfo.addRP(rpname, type=rptype, value=rpvalue, const=rpconst, range=rprange,
                             location=unitname, comment=rpcomment)
            # now for all documentation without corresponding parameter declarations
            # which start with __
            for (key,value) in list(unit['D'].items()):
                if key[:2] == '__' and key not in unit['PARAMETER']:
                   # a type of DOC is for documentation only
                   rpInfo.addRP(key,type='DOC',location=unitname, comment=value)
        # Generate enough plot_var_N names for all variables in the simulation.
        if type(max_plot_vars) == type(1):
            numRegularPlotVars = max_plot_vars
            for i in range(1,numRegularPlotVars+1):
                rpInfo.addRP('plot_var_' + '%d'%i, type="STRING", value='"none"',
                             location="IO/IOMain", comment="(automatically generated by setup)")
        else:
            rpInfo.addRP('plot_var_' + '%s'%max_plot_vars, type="STRING", value='"none"',
                         location="IO/IOMain", comment="(automatically generated by setup)")
        return rpInfo

    def getVarInfo(self,**kw):
        varInfo = VarInfo(**kw)
        for (unitname,unit) in list(self.units.items()):
            # Add all variables
            for (var,(attr,eos1,eos2)) in list(unit['VARIABLE'].items()):
                var = var.lower()
                varcomment = unit['D'].get("%s_variable" % var,"")
                varInfo.addVar(var,type="Variable", location=unitname,
                                   attribs=attr,comment=varcomment)
            # Add all particle properties
            for (var,attr) in list(unit['PARTICLEPROP'].items()):
                varcomment = unit['D'].get("%s_particleprop"%var,"")
                varInfo.addVar(var,type="Particle Property", location=unitname,
                                   attribs=[attr],comment=varcomment)
            # Add all particle types
            for var in list(unit['PARTICLETYPE'].keys()):
                varcomment = unit['D'].get("%s_particletype"%var,"")
                varInfo.addVar(var,type="Particle Type",location=unitname,
                                   attribs=[],comment=varcomment)
            # Add all particle property -> grid var maps
            for (var,attr) in list(unit['PARTICLEMAP'].items()):
                varInfo.addVar(var,type="Particle Property Map",
                                   location=unitname, attribs=[attr])
            # Add all species
            for var in list(unit['SPECIES'].keys()):
                varcomment = unit['D'].get("%s_species"%var,"")
                varInfo.addVar(var,type="Species",location=unitname,
                                   attribs=[],comment=varcomment)
            # Add all face variables
            for var in list(unit['FACEVAR'].keys()):
                varcomment = unit['D'].get("%s_facevar"%var,"")
                varInfo.addVar(var,type="Face Variable",location=unitname,
                                   attribs=[],comment=varcomment)
            # Add all scratch variables
            for var in list(unit['SCRATCHVAR'].keys()):
                varcomment = unit['D'].get("%s_scratchvar"%var,"")
                varInfo.addVar(var,type="Scratch Variable",location=unitname,
                                   attribs=[],comment=varcomment)
            for var in list(unit['SCRATCHCENTERVAR'].keys()):
                varcomment = unit['D'].get("%s_scratchcentervar"%var,"")
                varInfo.addVar(var,type="Scratch Center Variable",location=unitname,
                                   attribs=[],comment=varcomment)
            for var in list(unit['SCRATCHFACEXVAR'].keys()):
                varcomment = unit['D'].get("%s_scratchfacexvar"%var,"")
                varInfo.addVar(var,type="Scratch Face-X Variable",location=unitname,
                                   attribs=[],comment=varcomment)
            for var in list(unit['SCRATCHFACEYVAR'].keys()):
                varcomment = unit['D'].get("%s_scratchfaceyvar"%var,"")
                varInfo.addVar(var,type="Scratch Face-Y Variable",location=unitname,
                                   attribs=[],comment=varcomment)
            for var in list(unit['SCRATCHFACEZVAR'].keys()):
                varcomment = unit['D'].get("%s_scratchfacezvar"%var,"")
                varInfo.addVar(var,type="Scratch Face-Z Variable",location=unitname,
                                   attribs=[],comment=varcomment)
                
            # Add all Mass Scalars
            for var in list(unit['MASS_SCALAR'].keys()):
                varcomment = unit['D'].get("%s_mass_scalar"%var,"")
                varInfo.addVar(var,type="Mass Scalar",location=unitname,
                                   attribs=[],comment=varcomment)
            # Add all Fluxes
            for var in list(unit['FLUX'].keys()):
                varcomment = unit['D'].get("%s_flux"%var,"")
                varInfo.addVar(var,type="Flux",location=unitname,
                                   attribs=[],comment=varcomment)
        return varInfo


    # set default values for options possibly based of list of units we are using
    def adjustOpts(self):
        GVars.out.put("Computing default values for options not specified on command line",globals.IMPINFO)

        # Find all Grid related units we have now
        gP = globals.gridPrefix
        if not gP.endswith(os.sep): gP = gP+os.sep
        gridList = [ x[len(gP):] for x in self.getList() if x.startswith(gP) and x[len(gP):] ]
        # now gridList also includes non Grid subunits of "Grid/GridMain"
        # eliminate choices which do not occur in gridChoices
        grids = {}
        for x in gridList: # mark grids we have
            fp = x.split(os.sep)[0]
            if fp in globals.gridChoices: grids[fp] = 1
        if not grids:
           GVars.out.put("No Grid specified/found",globals.WARN)
           grid = ""
        elif len(grids) == 1:
           grid = list(grids.keys())[0]
           GVars.out.put("Using Grid %s" % grid,globals.INFO)
        else:
           raise SetupError("Found multiple grids %s!" % str(list(grids.keys())))

        # Adjust maxblocks
        if GVars.maxblocks==None: # using UG
           if grid in ('UG', 'UGSpatialFirst'):
             GVars.maxblocks = 1
           elif grid == 'AMR':  # using paramesh or Amrex or something like that
             if GVars.dimension==3: 	 
                GVars.maxblocks = 200
             else: 	 
                GVars.maxblocks = 1000
           else:    # some other grid package 
             if GVars.dimension==3: 	 
                GVars.maxblocks = 200 	 
             else: 	 
                GVars.maxblocks = 1000 	 

        # check if we are trying nonfbs and not UG
        if not GVars.setupVars.get("fixedBlockSize") and (grid not in ('UG', 'UGSpatialFirst', 'Chombo')):
           GVars.out.put("Non-Fixed Block size works only with UG",globals.ERROR)
           raise SetupError("Change your grid or switch to fixed block size")

        # check if we have included any particle units.
        GVars.particlesUnitIncluded = False
        for unit in self.getList():
            #We only check for Particles subdirectories because the top level
            #Particles directory (full of stubs) is included in ALL simulations.
            if unit.startswith("Particles/"):
                GVars.particlesUnitIncluded = True
                break

    def writeSetupUnitsFile(self):
        file = open(os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.SetupUnitsFilename), 'w')
        file.write("\n".join(self.getList()))
        file.write("\n")
        file.close()

    def generateSetupFlashUnits(self):
        fname = os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.SetupFlashUnitsFilename)
        tname = os.path.join(GVars.binDir,globals.SetupFlashUnitsTemplate)

        tpl = Template(tname)
        tpl["unit_names"] = self.getList()
        if len(tpl["unit_names"]) > 39:
            pUnits = self.getParentUnits()
            cUnits = [x for x in list(self.units.keys()) if x not in pUnits]
            if len(cUnits) <= 39: # if compression actually got us into the allowed range...
                cUnits.sort()
                GVars.out.put("Making list of units for %s more concise by removing parent units:\n   reduced %d -> %d unit names." \
                              % (fname,len(tpl["unit_names"]),len(cUnits)), globals.INFO)
                tpl["unit_names"] = cUnits
            else:
                 # exactly one Capitalized component after zero or more small components?
                weedable = re.compile("^([a-z][^/]*/)*([A-Z][^/]*)$")
                alwayskeep = re.compile("Main$") # never weed out a final *Main component
                cUnits = [x for x in cUnits if not weedable.match(x) or alwayskeep.match(x)]
                if len(cUnits) <= 39: # if more compression now got us into the allowed range...
                    cUnits.sort()
                    GVars.out.put("Making list of units for %s more concise by removing parent and stub-only units:\n   reduced %d -> %d unit names." \
                                  % (fname,len(tpl["unit_names"]),len(cUnits)), globals.INFO)
                    tpl["unit_names"] = cUnits
                else:
                    # exactly one Capitalized component after zero or more small components followed by "/localAPI"?
                    weedable = re.compile("^([a-z][^/]*/)*([A-Z][^/]*)/localAPI$")
                    alwayskeep = re.compile("Main/localAPI$") # never weed out a final *Main component
                    cUnits = [x for x in cUnits if not weedable.match(x) or alwayskeep.match(x)]
                    if len(cUnits) <= 39: # if more compression now got us into the allowed range...
                        cUnits.sort()
                        GVars.out.put("Making list of units for %s more concise by removing localAPI directories:\n   reduced %d -> %d unit names." \
                                      % (fname,len(tpl["unit_names"]),len(cUnits)), globals.INFO)
                        tpl["unit_names"] = cUnits
        tpl.generate(fname)

    def createMakefiles(self):
      """ Makefiles: one for each toplevel Unit. Makefiles for subunits get
      appended to the top one."""

      tname = os.path.join(GVars.binDir,globals.MakefileStubTemplate)
      tpl = Template(tname)
      for unitname in self.getList():
        lowestBase = getLowestBase(unitname)
        # if we find a unit (ie leaf is a capital
        # otherwise it's just an organizational directory
        # the exception is "flashUtilities" which is not a
        # Unit.  Look for it specifically
        if not is_upper(lowestBase[0]):
            if "flashUtilities" in lowestBase:
                lowestBase = "flashUtilities"
            else:
                continue

        source = os.path.join(GVars.sourceDir, unitname, "Makefile")
        target = 'Makefile.'+lowestBase
            
        if os.path.isfile(source):
            file = open(target, 'a')
            file.write(open(source).read())
            file.close()
        elif unitname==lowestBase:
            tpl["unitname"] = unitname
            tpl.generate(target)
            
################################################# Internally called code
# Return the part of name closest to ROOT dir, which starts with capital letter
def getLowestBase(base):
    parts = os.path.normpath(base).split(os.sep)
    x = [p for p in parts if is_upper(p[0])]
    try:
       return x[0]
    except IndexError:
       return base

def getDefs(dirName):
    defsList = []
    for f in os.listdir(dirName):
        fpath = os.path.join(dirName,f)
        if os.path.isfile(fpath) and os.path.splitext(fpath)[-1]==".ini":
            defsList.append(fpath )
    return defsList
