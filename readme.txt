# 2022/11/17
the TDI interp produce spikes, the shift 41 is to remove spikes in start transient, need to zero the spikes in end transient.


# 2022/08
this directory test spikes in XYZ (XYZ shift by 615 sec) and it's influence on snr and loglike. 
#
# there is no spike from td and fd shift, but 20220506 wh1 bestfit do have spike from shift.
# why I can't reproduce it?   check again.  20220820 
# (the spikes are from start/end transient resulting from TDI interp. The spikes in the end transient are delayed from plunge by TDI delay.)
#
When the ode evolve close to plunge, the td and fd  shift produce spikes in start and end transient.

