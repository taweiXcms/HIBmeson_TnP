This is ole summer 2014 version
TnP_HIBmeson
============

Tag and Probe for HI B meson

Usage:

1.
Run "reduce.cc" on Bfinder ntuple to reduce the file size.

2.
Run "TnPeff.cc" to create TnP efficiency result using side-band method and a TTree for later RooFit method use.

3.
Run "doTnPRooFit.cc", input = output of "TnPeff.cc", to create efficiency result using RooFit method.

4.
Run "DrawSF.cc", input = either output of "TnPeff.cc" or "doTnPRooFit.cc", to create the SF histograms.
