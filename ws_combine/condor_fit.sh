#!/bin/sh
i=$1 #Get array index from command line

ulimit -s unlimited
set -e
cd /afs/cern.ch/work/a/atarabin/SLtest/CMSSW_14_1_0_pre4/src
export SCRAM_ARCH=el9_amd64_gcc12
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd /afs/cern.ch/work/a/atarabin/SLtest/CMSSW_14_1_0_pre4/src/SLtest/ws_combine

mkdir fit_${i}
cd fit_${i}

echo "Running iteration $i"
new_input_file="../../inputs/htt_tt_125_8TeV_${i}.txt"
sed "s/htt_tt\.input\.root/htt_tt.input_${i}.root/g" ../../inputs/htt_tt_125_8TeV.txt > "$new_input_file"
text2workspace.py "$new_input_file" -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose --PO 'map=htt_tt_0_8TeV/ggH:r_high[1,-10,10]' --PO 'map=htt_tt_1_8TeV/ggH:r_low[1,-10,10]' --PO 'map=htt_tt_2_8TeV/ggH:r_low[1,-10,10]' -o ../../ws/ws_mu_ggH_${i}.root -m 125
combine -M MultiDimFit --algo grid -n .high.${i} -d ../../ws/ws_mu_ggH_${i}.root -P r_high
combine -M MultiDimFit --algo grid -n .low.${i} -d ../../ws/ws_mu_ggH_${i}.root -P r_low
plot1DScan.py higgsCombine.high.${i}.MultiDimFit.mH120.root --POI r_high -o scan_ggH_high_${i}
plot1DScan.py higgsCombine.low.${i}.MultiDimFit.mH120.root --POI r_low -o scan_ggH_low_${i}


