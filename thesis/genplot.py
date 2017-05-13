#!/usr/bin/python

# 2016 )c( by Tosten Lehmann

# This is a little script demonstrating how to plot to a pdf file
# using matplotlib (or pylab).
# Shows a few of common things to do with plots to get you started.
# With matplotlib, the possibilities are (almost) endless...

# Import all plotting and common matlab functionality
# (purists avert your eyes!).
from pylab import *

# Make sure to use type 1 (postscript) fonts
rcParams['ps.useafm']=True
rcParams['pdf.use14corefonts']=True
rcParams['font.family']='serif'
rcParams['font.serif']='Times'
rcParams['font.size']=10.0 # (pt) need to align with the figure size
rcParams['text.usetex']=True # needed for '-' sign

# Generate some data
seed(0)
N=200
timev=linspace(0,1,N)
sigv=sin(4*pi*timev)+(rand(N)-0.5)*0.5

# I find this size (5.5" by 3.5") good for including in documents
fig=figure(figsize=(5.5,3.5))

# In case you want the figure to go closer to the page margin than default
# and to adjust white space between subplots depending on label needs
subplots_adjust(left=0.11,bottom=0.11,right=0.89,top=0.97,wspace=0.3,hspace=0.2)
 
# Some common plotting functions
# (matlab-like commands should be self-explanatory)

ax1=subplot(111)
xlabel('Time [s]') # need to appear before twinx()
grid() # need to appear before twinx() for y-grid aligned with left axes
p1=plot(timev,sigv**2,color='#c5281c') #unsw red 2
ylabel('Power [W]')
ax2=ax1.twinx() # add another y-axes on the right
p2=plot(timev,sigv,color='#0092c8') #unsw blue 2
ylabel('Signal [V]')
legend((p1[0],p2[0]),('Power','Signal'))

# File name ending determines type - I like to use pdf figures
savefig('plotout.pdf')

