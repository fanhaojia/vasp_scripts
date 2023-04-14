#!/usr/bin/bash

git init 
git remote remove origin
git remote add origin git@github.com:fanhaojia/vasp_scripts.git
git init 
git add *
git commit -m "plot_soc_band_spin_projected"
git push origin master 


