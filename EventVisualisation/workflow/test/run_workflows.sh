#!/bin/bash
o2-eve-workflow --infile tpcdigits.root --input-type digits --output-type digits --monitoring-backend flume://127.0.0.1:7654 >& reco.log

