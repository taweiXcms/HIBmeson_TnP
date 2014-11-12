###Condor submitting template for plain root jobs. Run on all the files separately in a given folder, the DATASET folder below
###TA-WEI WANG, 02/20/2014 created
###Your plain root file needs to be modified , see example loop.C
###Please compile the .C before submission to make sure your code is working.
###Checking condor jobs status: condor_q <username> 

###Plain root .C to be run
CONFIGFILE="fitJpsi.C"

###All the header/related files needed
TRANSFERFILE="fitJpsi.C,utilities.h,format.h"

###Folder location within which files are to be run
#DATASET=/mnt/hadoop/cms/store/user/tawei/TnPBntuple/TnPnt_BoostedMC_20141022_Kp/*
DATASET=/mnt/hadoop/cms/store/user/tawei/TnPBntuple/TnPnt_BfinderBoostedMC_20141022_hckim-HIJINGemb_inclBtoPsiMuMu/*
#DATASET=/mnt/hadoop/cms/store/user/tawei/TnPBntuple/TnPnt_20141022_PAMuon_HIRun2013_28Sep2013_v1/*
#DATASET=/mnt/hadoop/cms/store/user/tawei/TnPBntuple/TnPnt_20141022_PAMuon_HIRun2013_PromptReco_v1/*

###Output file location
#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_BoostedMC_20140930_Kp
#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_TnPnt_BfinderBoostedMC_20141022_hckim-HIJINGemb_inclBtoPsiMuMu
#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_20140930_PAMuon_HIRun2013_28Sep2013_v1
#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_20140930_PAMuon_HIRun2013_PromptReco_v1

#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_TnPnt_BfinderBoostedMC_20141022_hckim-HIJINGemb_inclBtoPsiMuMu
DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_TnPnt_BfinderBoostedMC_20141022_hckim-HIJINGemb_inclBtoPsiMuMu_testAbsDxyDz
#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_TnPnt_BfinderBoostedMC_20141022_hckim-HIJINGemb_inclBtoPsiMuMu_testAbsDxyDz_test1
#DESTINATION=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/1022_DefaultMuId_TnPfiles/foutput_TnPnt_BfinderBoostedMC_20141022_hckim-HIJINGemb_inclBtoPsiMuMu_testAbsDxyDz_test2

###Output file name
OUTFILE="foutput"

###Maximum number of files to be run
MAXFILES=2000

###Log file location and it's name
LOGDIR=/net/hisrv0001/home/tawei/HeavyFlavor_20131030/TnP/HIBmeson_TnP/Code/logout
LOGNAME=testrootcondor

########################## Create subfile ###############################
dateTime=$(date +%Y%m%d%H%M)
fileCounter=0
INFILE=""
mkdir -p $DESTINATION
mkdir -p $LOGDIR

for file in $DATASET
do
  if [ $fileCounter -ge $MAXFILES ]
  then
  break
  fi
    INFILE="$file"
    fileCounter=$((fileCounter+1))
  if [ -f $DESTINATION/${OUTFILE}_${fileCounter}.root ]; then
    echo "Output already exists : ${OUTFILE}_${fileCounter}.root"
  else 
    # make the condor file
    cat > subfile <<EOF
    
Universe = vanilla
Initialdir = .
Executable = exec_condor.sh
+AccountingGroup = "group_cmshi.$(whoami)"
Arguments =  $CONFIGFILE $DESTINATION ${OUTFILE}_${fileCounter}.root $INFILE
GetEnv       = True
Input = /dev/null

# log files
Output       = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.out
Error        = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.err
Log          = $LOGDIR/$LOGNAME-$dateTime-${fileCounter}.log

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = $TRANSFERFILE

Queue
EOF
############################ Submit ###############################
    #cat subfile
    condor_submit subfile
    mv subfile $LOGDIR/$LOGNAME-$dateTime-$fileCounter.subfile
  fi
done
echo "Submitted $fileCounter jobs to Condor."
