#!/bin/bash
#
# *** SUBMIT ENSEMBLE ***
#

# specify base-config filename
BASECONFIG='myBASECONFIG'
# specify path to any user-configs ('/' otherwise)
USERCONFIGPATH='myDIR'
# specify experiment (user-config) filename, excluding the .xy ensemble member extension
ENSEMBLEID='myENSEMBLEMEMBER'
# specify run duration (integer years)
YEARS='10000'
# specify any restart name (empty string otherwise)
RESTARTID=''
# set first parameter axis max # (0-9)
memberimax=9
# set second parameter axis max # (0-9)
memberjmax=9
# specific any particular queue to be used (empty string otherwise)
QUEUE='-q cat.q'

# initialize loop counter i
memberi=0
while [ $memberi -le $memberimax ]; do
  # initialize loop counter j
  memberj=0
  while [ $memberj -le $memberjmax ]; do

    # set index and userconfig name
    RESTART=$RESTARTID
    EXPERIMENT=$ENSEMBLEID"."$memberi$memberj
  
    # submit
    echo $EXPERIMENT"/"$RESTART
    qsub $QUEUE -j y -o cgenie_log -V -S /bin/bash runmuffin.sh $BASECONFIG $USERCONFIGPATH $EXPERIMENT $YEARS
    sleep 10
    #qstat -f

    #
    let memberj=$memberj+1
  done
  let memberi=$memberi+1
done

#
qstat -f
    