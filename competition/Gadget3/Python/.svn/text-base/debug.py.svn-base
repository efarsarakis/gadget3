# --------------------------
#           USAGE
#---------------------------

# 1. attach gdb to running program (e.g.: using the debug_run() function and the 
#    PAUSE_RUN_TO_ATTACH_DEBUGGER Config.sh option in Gadget) or use it to
#    load a core dump (e.g.: for a single file by "gdb P-Gadget3 core" or for multiple files
#    by giving the core_base variable to gravtree() or particle_data())
#    you need a fairly new gdb version, the code was tested with gdb version 7.4.50
# 
# 2a. type in gdb: "py execfile("Python/debug.py")" this loads this file
#     then load, e.g., the tree by typing " py tree = gravtree() " within gdb
#     the tree is then loaded and can be analysed by standard Python commands
#     note that when typing them within gdb you always have to type "py" at the beginning of the line
#
# 2b. or to work more comfortably with ipython to the following instead of the instructions under (2a.)
#     (you need a new version of ipython >= 0.12, compiled with support for pyzmq)
#     type in gdb: " py execfile("Python/gdb_ipython.py") "
#     this will output a line like " [IPKernelApp] --existing kernel-19836.json "
#     then open another console and start " ipython console --existing kernel-19836.json "
#     of course using the kernel number that was printed in the gdb window
#     in ipython then type " execfile("Python/debug.py") " and continue as under point (2a.)
#     in my test connceting IPython to this kernel was, however, somewhat unstable
#
# currently only tested when turning off optimization (-O0 -g) when compiling Gadget 

#---------------------------

import gdb

from pylab import *

#---------------------------

# USING THE debug_run CALSS FOR INTERACTIVE PARALLEL DEBUGGING
#
# 1. compile Gadget with PAUSE_RUN_TO_ATTACH_DEBUGGER option
# 2. start your Gadget MPI run (currently all tasks need to be
#    on the same node
# 3. start "gdb P-Gadget3"
# 4. type in gdb: "py execfile("Python/debug.py")"
# 5. type in gdb: "myrun = debug_run()" -> gdb attached to MPI tasks
# 6. set, e.g., a breakpoint by "break forcetree.c:2449"
# 7. start run with "py myrun.start()"
# 8. once the breakpoint is reached by all task you can start debugging
# 9. change between MPI tasks by typing e.g. "inferior 3", which
#    changes to task 3
# 10. or load data from all tasks e.g. with "py part = particle_data()"
# 11. continue run with "py myrun.cont()"

class debug_run:

  def __init__ (self):  # attaches to all tasks of paused Gadget run (need to be on same node) 
 
    self.pidlist = loadtxt("pid_list_for_debugger.txt",dtype=int32)
    print "PID list loaded:"
    print self.pidlist
  
    self.NTask = self.pidlist.size
    print "NTask = ", self.NTask
 
    gdb.execute("set target-async 1")
    gdb.execute("set pagination off")
    gdb.execute("set non-stop on")   # manage tasks individually
 
    if self.NTask>1:
      gdb.execute("clone-inferior -copies " + str(self.NTask-1))
    
    for task in range(self.NTask):
      gdb.execute("inferior " + str(task+1))
      gdb.execute("attach " + str(self.pidlist[task]))
      gdb.execute("interrupt")
      
    #gdb.execute("set schedule-multiple")  # to manage all tasks at the same time
      
  def start(self):  # starts running
  
    for task in range(self.NTask):
      gdb.execute("inferior " + str(task+1))
      
      str_backtrace = (gdb.execute("backtrace", to_string=True)).splitlines()
      for i in range(str_backtrace.__len__()):   # steps out of sleep function
        if str_backtrace[i].find("pause_run_to_attach_debugger") > 0:
          for j in range(i):
            gdb.execute("up")
          break
            
      gdb.execute("set var continue_run = 1")
      
    gdb.execute("continue -a")
    
  def cont(self):  # continue running
      
    gdb.execute("continue -a")
    
  def where(self):
    for task in range(self.NTask):
      gdb.execute("inferior " + str(task+1))
      gdb.execute("where")
      
  def inter(self):
    for task in range(self.NTask):
      gdb.execute("inferior " + str(task+1))
      gdb.execute("interrupt")
      
  def is_running(self):  # check if at least one task is still running
    
    for task in range(self.NTask):
      print task
      gdb.execute("inferior " + str(task+1))
      if gdb.selected_thread().is_running():
	return True
    
    return False

#---------------------------

def get_data_type (gdb_field_type):   # returns data type and if an array number of elements for struct field 
  
  fieldsize = gdb_field_type.sizeof   # find size in bytes of field
  typename = str(gdb_field_type)   # name of data type
  
  if typename.find("[") >= 0:   # is an array
    arrsize = int32(typename[typename.find("[")+1:typename.find("]")])
  else:
    arrsize = 1

  if typename.find("MyFloat")>=0:   # find data type for field
    if fieldsize/arrsize == 4:
      dattype = float32
    elif fieldsize/arrsize == 8:
      dattype = float64
    else:
      assert False
  
  elif typename.find("MyDouble")>=0:   # find data type for field
    if fieldsize/arrsize == 4:
      dattype = float32
    elif fieldsize/arrsize == 8:
      dattype = float64
    else:
      assert False
  
  elif typename.find("MyLongDouble")>=0:   # find data type for field
    if fieldsize/arrsize == 4:
      dattype = float32
    elif fieldsize/arrsize == 8:
      dattype = float64
    elif fieldsize/arrsize == 16:
      dattype = float128
    else:
      assert False
  
  elif typename.find("double")>=0:
    dattype = float64
  
  elif typename.find("float")>=0:
    dattype = float32
  
  elif typename.find("integertime")>=0:
    if fieldsize/arrsize == 4:
      dattype = int32
    elif fieldsize/arrsize == 8:
      dattype = int64
    else:
      assert False
  
  elif typename.find("short unsigned int")>=0: 
    dattype = uint16
  
  elif typename.find("unsigned int")>=0: 
    dattype = uint32
  
  elif typename.find("short int")>=0: 
    dattype = int16
  
  elif typename=="int" or typename.find("int ")>=0: 
    dattype = int32
  
  elif typename.find("MyIDType")>=0:
    if fieldsize/arrsize == 4:
      dattype = uint32
    elif fieldsize/arrsize == 8:
      dattype = uint64
    else:
      assert False
  
  else:
    print "type not known by get_data_type():", typename
    assert False

  return dattype, arrsize, fieldsize

#---------------------------

def get_field_data (item):   # get names, offsets, fieldsizes, and data types of all struct fields and fiels of nested structs & unions
  fields = item.type.fields()   # get fields
  fieldnum = fields.__len__()   # get number of fields
  
  names = []
  offsets = []
  dattypes = []
  fieldsizes = []
  arrsizes = []
  
  for i in range(fieldnum): 
    if str(fields[i].type).find("union")<0 and str(fields[i].type).find("struct")<0:   # normal field, i.e. no union or struct
      offset = uint64(item[fields[i].name].address)- uint64(item.address)   # find memory offset of struct field
      dattype,arrsize,fieldsize = get_data_type(fields[i].type)
      
      names.append(fields[i].name)
      offsets.append(offset)
      dattypes.append(dattype)
      fieldsizes.append(fieldsize)
      arrsizes.append(arrsize)
      
    else:   # if struct or union
      struct_offset = uint64(item[fields[i].name].address)- uint64(item.address)   # find memory offset of sub-struct
    
      name,offset,fieldsize,dattype,arrsize = get_field_data(item[fields[i].name])   # get content of sub-struct
      for j in range(name.__len__()):
        name[j] = fields[i].name + "." + name[j]
        offset[j] += struct_offset
      
      names.extend(name)   # add results of sub-struct
      offsets.extend(offset)
      dattypes.extend(dattype)
      fieldsizes.extend(fieldsize)
      arrsizes.extend(arrsize)

  return names, offsets, fieldsizes, dattypes, arrsizes

# ----- class for tree -----

class gravtree:
  
  def __init__ (self, core_base="none"):   # loads tree, give first part of core file names to loop over them
    print
    print "loading tree into Python"
    
    targetnum = 0
    core_names = []
    
    if (core_base=="none"):
      targetnum = gdb.inferiors().__len__()
    else:
      import os
      dirlist = os.listdir(".")
      for cur_file in dirlist:
        if cur_file.find(core_base)>=0: # found suitable core fike
          targetnum += 1
          core_names.append(cur_file)
    
    self.Numnodestree_tasks = zeros(targetnum,dtype=int32)
    self.MaxNodes_tasks = zeros(targetnum,dtype=int32)
    self.NTopnodes_tasks = zeros(targetnum,dtype=int32)
    self.NTopleaves_tasks = zeros(targetnum,dtype=int32)
    self.tasknums = zeros(targetnum,dtype=int32)
    self.task_inds = dict()
    
    for j in range(targetnum):
      if (core_base=="none"):
        print "doing task", j
        gdb.execute("inferior " + str(j+1)) # go to right task
      else:
        print "doing core file", core_names[j]
        gdb.execute("target core "+core_names[j]) # select right core file
    
      infer = gdb.selected_inferior()
    
      if j==0:
        self.MaxPart = uint64(gdb.parse_and_eval("All.MaxPart"))
    
      Numnodestree = uint64(gdb.parse_and_eval("Numnodestree"))
      MaxNodes = uint64(gdb.parse_and_eval("MaxNodes"))
      NTopnodes = uint64(gdb.parse_and_eval("NTopnodes"))
      NTopleaves = uint64(gdb.parse_and_eval("NTopleaves"))
      print "number of nodes = ", Numnodestree
      print "max number of nodes = ", MaxNodes
      print "number of topnodes = ", NTopnodes
      print "number of topleaves = ", NTopleaves
  
      self.task_inds[j] = arange(self.Numnodestree_tasks.sum(),self.Numnodestree_tasks.sum()+Numnodestree, dtype=uint32)   # index of current task in nodes data, for accessing data of specific task
      
      self.tasknums[j] = int32(gdb.parse_and_eval("ThisTask"))   # task numbers of tasks
      self.Numnodestree_tasks[j] = Numnodestree
      self.MaxNodes_tasks[j] = MaxNodes
      self.NTopnodes_tasks[j] = NTopnodes
      self.NTopleaves_tasks[j] = NTopleaves
  
      first_node = gdb.parse_and_eval("Nodes_base[0]")
      second_node = gdb.parse_and_eval("Nodes_base[1]")
      
      if j==0:
        self.nodesize = uint64(second_node.address) - uint64(first_node.address)
        print "sizeof(NODE) = ", self.nodesize
    
        names,offsets,fieldsizes,dattypes,arrsizes = get_field_data(first_node)   # get fields of tree node
        self.fieldnum = names.__len__()  # get number of fields 
        
        self.nodes = dict()   # dictionary for node data
      
        self.NTask = int32(gdb.parse_and_eval("NTask"))
   
      assert self.NTask == int32(gdb.parse_and_eval("NTask")) 
   
      rawdat = array(bytearray(infer.read_memory(first_node.address, Numnodestree * self.nodesize)))   # reads tree data
      rawdat = rawdat.reshape((Numnodestree, self.nodesize))   # reshapes the array so that the data of one node is in each row   
    
      for i in range(self.fieldnum):   # load main struct data, note that topnodes are loaded several times, i.e. for each task 
        dat = (rawdat[:,offsets[i]:offsets[i]+fieldsizes[i]].copy()).view(dtype=dattypes[i])   # extracts data for field 
        
        if j==0:
          self.nodes[names[i]] = dat
        else:
          self.nodes[names[i]] = append(self.nodes[names[i]], dat, axis=0)

      del rawdat

    self.get_bitflags()   # decodes bitflags

    self.keylist = self.nodes.keys()
    self.keylist.sort()
    print "fields read:"
    for i in range(self.keylist.__len__()): 
      print "  ",self.keylist[i]
    print

    self.targetnum = targetnum
    self.taskstart = append(array([0],dtype=int32), self.Numnodestree_tasks.cumsum()[:-1]) # top node of each task
    
    if self.NTask != targetnum:
      print "WARNING not all core files loaded!"
      print "NTask =", self.NTask, ", number of core files =", targetnum

  def get_bitflags(self):   # sets flags of node
    self.nodes['TOPLEVEL'] = array((self.nodes['u.d.bitflags'] >> 0) % 2, dtype=bool)
    self.nodes['DEPENDS_ON_LOCAL_MASS'] = array((self.nodes['u.d.bitflags'] >> 1) % 2, dtype=bool)  
    self.nodes['MAX_SOFTENING_TYP'] = (self.nodes['u.d.bitflags'] >> 2) % 8
    self.nodes['MIXED_SOFTENINGS_IN_NODE'] = array((self.nodes['u.d.bitflags'] >> 5) % 2, dtype=bool)
    self.nodes['INTERNAL_TOPLEVEL'] = array((self.nodes['u.d.bitflags'] >> 6) % 2, dtype=bool)
    self.nodes['MULTIPLEPARTICLES'] = array((self.nodes['u.d.bitflags'] >> 7) % 2, dtype=bool)
    self.nodes['NODEHASBEENKICKED'] = array((self.nodes['u.d.bitflags'] >> 8) % 2, dtype=bool)
    self.nodes['INSIDE_LINKINGLENGTH'] = array((self.nodes['u.d.bitflags'] >> 9) % 2, dtype=bool)
    
  def get_tasks(self):
    self.nodes['task'] = -1*ones(self.Numnodestree_tasks.sum(), dtype=int32)
    self.nodes['taskind'] = -1*ones(self.Numnodestree_tasks.sum(), dtype=int32)
    
    for i in range(self.NTask):
      self.nodes['task'][self.task_inds[i]] = self.tasknums[i]
      self.nodes['taskind'][self.task_inds[i]] = i
    
  def get_node_level(self, node):
    taskind = self.nodes['taskind'][node]
    if node==self.taskstart[taskind]:
      return 1
    elif self.nodes['levels'][node] >= 0:
      return self.nodes['levels'][node]
    else:
      return self.get_node_level(self.nodes['u.d.father'][node] - self.MaxPart + self.taskstart[taskind]) + 1
    
  def get_levels(self):   # find tree level for each node, 0 corresponds to top node covering the whole simulation volume 
    print "getting node levels ..."
    
    try:
      self.nodes['task']
    except:
      self.get_tasks()
    
    self.nodes['levels'] = -1*ones(self.Numnodestree_tasks.sum(), dtype=int32)
    for i in range(self.Numnodestree_tasks.sum()):
      self.nodes['levels'][i] = self.get_node_level(i)
    print "done!"
    self.maxlevel = self.nodes['levels'].max()
    print "maximum level of tree =", self.maxlevel 
    
  def show_stats(self):   # show some stats of tree
    try:
      self.maxlevel
    except:
      self.get_levels()
      
    for i in arange(self.maxlevel+1)[1:]:
      ind  = where(self.nodes['levels']==i)[0]
      num_topnodes = where(self.nodes['TOPLEVEL'][ind])[0].size/self.targetnum
      print "level", i, " nodes", ind.size - num_topnodes*(self.targetnum-1), " topnodes", num_topnodes

  def get_index_unique(self):   # excludes duplicated topnodes
    for i in range(self.NTask):
      if i == 0:
        self.ind_unique = self.task_inds[0]
      else:
        self.ind_unique = append(self.ind_unique, self.task_inds[i][self.NTopnodes_tasks[i]:])

# ----- class for particles -----

class particle_data:
  def __init__ (self, core_base="none"): # loads particles, give first part of core file names to loop over them
    
    print
    print "loading particles into Python"
    
    targetnum = 0
    core_names = []
    
    if (core_base=="none"):
      targetnum = gdb.inferiors().__len__()
    else:
      import os
      dirlist = os.listdir(".")
      for cur_file in dirlist:
        if cur_file.find(core_base)>=0: # found suitable core fike
          targetnum += 1
          core_names.append(cur_file)
    
    self.NumPart_tasks = zeros(targetnum,dtype=int32)
    self.N_gas_tasks = zeros(targetnum,dtype=int32)
    self.tasknums = zeros(targetnum,dtype=int32)
    self.task_inds = dict()
    
    for j in range(targetnum):
      if (core_base=="none"):
        print "doing task", j
        gdb.execute("inferior " + str(j+1)) # go to right task
      else:
        print "doing core file", core_names[j]
        gdb.execute("target core "+core_names[j]) # select right core file
    
      infer = gdb.selected_inferior()
    
      NumPart = uint64(gdb.parse_and_eval("NumPart"))
      N_gas = uint64(gdb.parse_and_eval("N_gas"))
      MaxPart = uint64(gdb.parse_and_eval("All.MaxPart"))
      print "number of particles = ", NumPart
      print "number of gas particles = ", N_gas
      print "max number of particles = ", MaxPart
      
      self.task_inds[j] = arange(self.NumPart_tasks.sum(),self.NumPart_tasks.sum()+NumPart)   # index of current task in part data, for accessing data of specific task
      
      self.tasknums[j] = int32(gdb.parse_and_eval("ThisTask"))   # task numbers of tasks
      self.NumPart_tasks[j] = NumPart
      self.N_gas_tasks[j] = N_gas
      
      if j==0:
        self.MaxPart = MaxPart
      else:
        assert self.MaxPart == MaxPart
      
      first_part = gdb.parse_and_eval("P[0]")
      second_part = gdb.parse_and_eval("P[1]")
      
      if j==0:
        self.partsize = uint64(second_part.address) - uint64(first_part.address)
        print "sizeof(PARTICLE) = ", self.partsize

        names,offsets,fieldsizes,dattypes,arrsizes = get_field_data(first_part)
        self.fieldnum = names.__len__()
        
        self.part = dict()   # dictionary for node data
        
        self.NTask = int32(gdb.parse_and_eval("NTask"))
        
      assert self.NTask == int32(gdb.parse_and_eval("NTask"))

      rawdat = array(bytearray(infer.read_memory(first_part.address, NumPart * self.partsize)))   # reads particle data
      rawdat = rawdat.reshape((NumPart, self.partsize))   # reshapes the array so that the data of one particle is in each row   
      
      for i in range(self.fieldnum):   # load main struct data
        dat = (rawdat[:,offsets[i]:offsets[i]+fieldsizes[i]].copy()).view(dtype=dattypes[i])   # extracts data for field 
        if j==0:
          self.part[names[i]] = dat
        else:
          self.part[names[i]] = append(self.part[names[i]], dat, axis=0)
      
      del rawdat
    
    self.keylist = self.part.keys()
    self.keylist.sort()
    print "fields read:"
    for i in range(self.keylist.__len__()): 
      print "  ",self.keylist[i]
    print
    
    self.targetnum = targetnum
    
    if self.NTask != targetnum:
      print "WARNING not all core files loaded!"
      print "NTask =", self.NTask, ", number of core files =", targetnum
