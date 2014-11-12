CONFIGFILE=$1
DESTINATION=$2
OUTFILE=$3
INFILE=$4

source /osg/app/cmssoft/cms/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc472
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc472/cms/cmssw/CMSSW_5_3_20/
eval `scramv1 runtime -sh` 
cd -

root -l -b -q  $CONFIGFILE++\(\"${INFILE}\",\"${OUTFILE}\"\)
#root -l -b -q  $CONFIGFILE\(\"${INFILE}\",\"${OUTFILE}\"\)
#root -l -b -q  $CONFIGFILE++
mv ${OUTFILE} ${DESTINATION}/${OUTFILE}
