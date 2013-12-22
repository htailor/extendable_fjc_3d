Extendable Freely Jointed Chain 3D
==================================

Model of the extendable freely-jointed-chain in three dimensions. 

Using the Fourier Trasform Integral method, the script plots the Partition Function for an Extendable Freely Jointed Chain Model at different extensions. For a freely-jointed-chain with `N` links, each of length `a`, we include a parameter `p` that controls the extensibility of each link.

This calculation supports the Extendable Freely Jointed Chain chapter of the UCL PhD thesis.

Requirements
============

Python 2.7
TexLive distribution

Has been tested on Mac OS X, Windows 8 and Linux

Instructions
============

When the script is executed a `results` directory is created which contains all the raw data and plots (pdf format).

To run different parameters simple change the variables `a` and `p` near the top of the script. In the main section you can specify the chain size in the `N` list. 

**WARNING: Using N > 12 causes unstable results that will increase run time (see thesis for more info)**


Files
=====

`efjc_3d.py` - Main file to execute

`libmath.py`
