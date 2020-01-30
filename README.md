## Bruteforce attack tool for a Legendre based prng
This repository contains a set of python CLI tools to launch a bruteforce attack towards a Legendre symbol based PRF, 
as the one for which Ethereum offers [https://legendreprf.org/bountyinstances](bounties), alongside with a script that 
produces charts from data collected by several bruteforce instances.
- pure brutefore can be done via `legendre_prng.py`, usage info with `python3 legendre_prng.py --help 
- charts can be produced (and found inside the `figs` folder) by launching `python charts.py`. Default parameters can 
cause the execution to take quite long (~10 mins) and can be changed in the file itself.
