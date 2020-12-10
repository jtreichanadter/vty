#!/usr/bin/env bash
# ---------------------------------------------------------------------------- #
# Adjust command-line environment for 'vty' tools
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #



# load required modules (potentially unnecessary)
echo "checking modules"
if command -v module; then
  module load python  ||  echo -e "\033[33mCommand 'module load python' threw an error.  Possibly okay...\033[0m"
  module load gnuplot ||  echo -e "\033[33mCommand 'module load gnuplot' threw an error. Possibly okay...\033[0m"
fi


# show utility versions
if command -v python; python --version
if command -v gnuplot; gnuplot -V


# add vty library to path
PATH="~/.scripts/vty:${PATH}"
export PATH


# set python environment alias
alias pyvty="python -B -i ~/.scripts/vty/vty.py"






# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# author    = (Jonathan) Tyler Reichanadter
# email     = jtreichanadter@berkeley.edu
# copyright = Copyright 2020, Neaton Group at UC Berkeley
# license   = GPL
