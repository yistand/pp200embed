#! /bin/tcsh

#starver SL11d
starver SL12d

cons StMiniMcTree+

set nevents="-1"
root4star -b -q minimc_macro.cc\($nevents,\"lists/ppProton200.txt\",\"Proton\",200\)
root4star -b -q minimc_macro.cc\($nevents,\"lists/ppAntiProton200.txt\",\"AntiProton\",200\)
root4star -b -q minimc_macro.cc\($nevents,\"lists/ppKplus200.txt\",\"Kplus\",200\)
root4star -b -q minimc_macro.cc\($nevents,\"lists/ppKminus200.txt\",\"Kminus\",200\)
root4star -b -q minimc_macro.cc\($nevents,\"lists/ppPiplus200.txt\",\"Piplus\",200\)
root4star -b -q minimc_macro.cc\($nevents,\"lists/ppPiminus200.txt\",\"Piminus\",200\)

