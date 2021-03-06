#!/usr/bin/env Rscript

cd "$(dirname "$0")"
PC=$(hostname) # add user name to log for avoiding clashes between machines
LOG_FILE="${PC}_simulation.log"
Rscript simulation.R | tee ${LOG_FILE}
