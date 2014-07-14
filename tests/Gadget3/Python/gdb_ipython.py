from pylab import *

from IPython.zmq.ipkernel import IPKernelApp
 
app = IPKernelApp.instance()
app.initialize([])
app.start()
