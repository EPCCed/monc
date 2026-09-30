#!/bin/bash

# ------------------------------------------------------------------------------------- #
#  Routine to check on all (or a subset) test_harness configurations.
#  Checks running status and stdout files for a few kinds of errors.
#  Optionally submits jobs that have not been submitted or resubmits those that have
#  encountered errors.
#
#  NOTE: This does NOT handle the "standard" cases.
#        Mainly because they all have the same job name and have more complex paths.
#        Ultimately, it might be simpler to use the existing python dictionaries to work
#        on those jobs (and perhaps these ones, too).
#
#  Input (all optional):
#    -a  Perform submission and restart actions.
#          DEFAULT: Not enabled. Only checks on status and logs success.
#    -c  Provide a single text pattern that will be globbed as:
#              test_harness/*${pattern}*
#        to select configuration directories for evaluation.
#          DEFAULT: All existing configuration directories are evaluated
#    -r  Number of restarts to allow for each job.
#          DEFAULT: 2
#    -o  Don't ask to confirm settings.
#          DEFAULT: Confirm settings with user.
#   
#
#  Output:
#    stdout - print to screen job status ['still running', 'ATTENTION', or 'Success']
#           - tail of most recent cycle stdout
#    copies of stdout from failed jobs as case.rerun in test_harness/.../monc_stdout
#    empty case.success files for completed cases as case.success in test_harness/.../monc_stdout
#
#  Actions:
#    Does nothing for running jobs
#    Marks successful jobs with a .success file in the monc_stdout directory
#      - .success files will contain any 'Miss match or un-enabled' warnings from the
#        final submission cycle
#    Attempts to restart failed jobs (two attpemts only, unless changed with -r)
#      - on second attempt, the previous checkpoint is deleted to restart further back
#      - If you wish to try more times, you can clear the contents of the offending .rerun file,
#        being careful to clear checkpoints as needed, too.
#
#  Example usage:
#    To run on all cases without submission actions:
#      ck_progress.sh
#    To run with submission and restart actions enabled for a very specific configuration
#    directory pattern:
#      ck_progress.sh -a -c ARCHER2_main_component_gnu_high
# ------------------------------------------------------------------------------------- #


RED='\033[0;31m'   # styles: 0:normal, 1:bold, 4:underline
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[0;32m'
UPurple='\033[4;35m'
NC='\033[0m' # No Color
owd=$(pwd)
cd $owd
if [ ! -d ${owd}/test_harness ] ; then
  echo "  STOP - test_harness directory not present"
  echo "  Execute from a MONC main directory"
  exit 1
fi

# Configure defaults
# We'll check on all test_harness runs beginning with A or M, following the naming
# conventions of the monc_test_harness.py 
list=$(/bin/ls -d test_harness/[A,M]*)

# Handle options
print_usage() {
  echo ""
  echo "-------------------------------------------------------------------------"
  echo "Usage of ck_progress.sh"
  echo "  Example usage: ck_progress.sh -a -c ARCHER2_main_component_gnu_high"
  echo "    - This usage runs with submission and restart actions enabled "
  echo "      and only one configuration will be evaluated."
  echo ""
  echo "Options:"
  echo "  -a  Perform submission and restart actions."
  echo "        DEFAULT: Not enabled. Only checks status and logs success."
  echo "  -c  Provide a single text pattern that will be globbed as:"
  echo "         test_harness/*\${pattern}*"
  echo "      to select configuration directories for evaluation."
  echo "        DEFAULT: All existing configuration directories are evaluated."
  echo "  -r  Number of restarts to allow for each job."
  echo "        DEFAULT: 2"
  echo "  -o  Don't ask to confirm script parameters."
  echo "        DEFAULT: Confirm script parameters with user."
  echo "-------------------------------------------------------------------------"
  echo ""
}

# Initialise restart limit
#trjh (because we subtract 1 later, this runs into trouble early, not deleting multiple checkpoints) rlimit=2
rlimit=3

while getopts 'ac:r:oh' flag; do
  case "${flag}" in
    a) act=1 ;;
    c) pattern="${OPTARG}" ;;
    r) change_rlimit="${OPTARG}" ;;
    o) override=1 ;;
    h) print_usage ; exit ;;
    *) print_usage ; exit 1 ;;
  esac
done

# Messages for options
if [[ ! -z $act ]] ; then
  echo -e "\n${RED}Performing submission & restart actions.${NC}\n"
else
  echo -e "\n${YELLOW}Running without submission & restart actions.${NC}"
  echo -e "${YELLOW}Action messages are still provided to show which actions would be taken.${NC}\n"
fi

if [[ ! -z $pattern ]] ; then
  echo " -- Overriding default directories -- "
  list=$(/bin/ls -d test_harness/*${pattern}*)
fi
echo -e "${BLUE}Evaluating:${NC}"
for inc in $list ; do
  echo "  $inc"
done

if [[ ! -z $change_rlimit ]] ; then
  echo -e "\n${GREEN} Changing number of restarts from 2 to: ${change_rlimit} ${NC}"
  rlimit=${change_rlimit}
fi

echo ''

# User verification of setup parameters
if [[ -z $override ]] ; then
  echo "Do you wish to proceed with these parameters?"
  select yn in "Yes" "No"; do
    case $yn in
      Yes ) break ;;
      No ) exit ;;
    esac
  done
fi

# Initialise number of jobs in queue for this user
niq=0
rtest=$(( $rlimit - 1 ))

# Configure submission and commands based on local machine scheduler
if [ -x "$(command -v qsub)" ] ; then
  system='pbs'
  submit="qsub"
  jn_search="PBS -N"
  jn_cut= | rev | cut -f 1 -d ' ' | rev
  qstat="/opt/ukmo/supported/bin/qstat_snapshot -u ${USER}"
  qstat_offset=5
  run_cond=" R \| Q "
  held_cond=" H  "
  jselect="qselect -u $USER -N"
  nqlim=200  # value for MONSOON normal queue
  job_deletion() {
    qselect -u $USER -N "$1" | xargs qdel
  }
elif [ -x "$(command -v sbatch)" ] ; then
  system='slurm'
  submit="sbatch"
  jn_search="SBATCH --job-name="
  jn_cut=" | cut -f 2 -d '='"
  qstat="squeue -l -u ${USER} --Format=jobid:10,state:10,name:80,reasonlist:15"
  qstat_offset=2
  run_cond="PENDING\|RUNNING\|CONFIGUR"
  held_cond="(Dependency)"
  jselect="sacct -n -X --format jobid --name"
  nqlim=128 # value for ARCHER2 in taskfarm QoS
  job_deletion() {
    scancel -n "$1"
  }
else
  echo "Error.  Unknown batch submission protocol."
  exit
fi



################################## Work through the list #######################################
for inc in $list ; do 
  echo "------- ------- ------ ------ ----- ----- ---- ---- --- --- -- -- - -  -  -   -   -    -"
  echo -e "${UPurple}${inc}${NC}"
  cd $inc/monc_stdout
  cwd=`pwd`  # We'll default to working from the local monc_stdout

  # Obtain case list
  cd ../tmp_config_dir
  cases=$(/bin/ls * | sort)
  cd ${cwd}

  # Iterate over cases
  for cnc in $cases ; do 

    # Obtain job name for this case
    if [[ "${system}" == "pbs" ]] ; then
      jn=$(grep "${jn_search}" ../tmp_qsub_dir/submonc_${cnc}  | rev | cut -f 1 -d ' ' | rev)
    elif [[ "${system}" == "slurm" ]] ; then
      jn=$(grep "${jn_search}" ../tmp_qsub_dir/submonc_${cnc} | cut -f 2 -d '=')
    fi
    echo -e "${BLUE}${jn}${NC}"

    # Check for already successful completion
    if [ -f ${jn}.success ] ; then
      echo -e "${GREEN} .success ${NC}"
      echo ''
      continue # check next case
    fi

    # Get most recent stdout for this case
    lfile=$(/bin/ls -v ou*$cnc* 2> /dev/null | tail -1) 

    # Get number of jobs in queue for this user
    # Disable actions if needed.
    niq=$(( $(eval $qstat | wc -l) - $qstat_offset ))
    if [[ ! -z $act ]] && [ ${niq} -ge ${nqlim} ] ; then
      echo -e "${YELLOW}\n    NUMBER IN QUEUE: ${niq} of ${nqlim} limit.${NC}"
      echo -e "${RED}    ---  DISABLING ACTIONS ---${NC}"
      act=
    fi


    # Check if still running, continue if so (don't act on running/queued jobs)
    if eval $qstat | grep -w ${jn} | grep -q "${run_cond}" ; then
      echo -e "${YELLOW}  ...still running...${NC}"
      if [ -z $lfile ] ; then
        echo '  There is no stdout for this case.'
      else
        echo ${cwd}/$lfile
        tail -3 $lfile
        grep -i 'model time' $lfile | tail -1
      fi
      echo ""
      continue # check next case
    fi

    # Handle situation where no output files for cases not queued/running
    if [ -z $lfile ] ; then
      echo '  There is no stdout for this case.'
      if [[ ! -z $act ]] ; then
        cd $owd
        if [ ${niq} -lt ${nqlim} ] ; then
          ${submit} ${inc}/tmp_qsub_dir/submonc_${cnc}
        else
          echo -e "${YELLOW}  -- Queue limit exceeded.  Job not submitted. --  ${NC}"
        fi
        cd $cwd 
      fi
      echo ""
      continue  # check next case
    fi

    # Print stdout tail for information
    echo ${cwd}/$lfile
    tail -3 $lfile 

    # Grab jobid information
    if grep -q 'This cycle job:' $lfile ; then
      thisjob=`grep 'This cycle job:' $lfile | cut -f 2 -d ':'`
      nextjob=`grep 'Next cycle job:' $lfile | cut -f 2 -d ':'`
    fi

    # For cases that don't have running jobs, inspect available information.

    # Check for errors
    if grep -qi 'error\|Caught signal\|ATP analysis\|Segmentation' $lfile ; then
      echo ''
      echo -e "${RED}  !!! found error !!!${NC}"
      echo ''
      echo $cwd/$lfile
      grep 'error\|Caught signal\|ATP analysis\|Segmentation' $lfile | uniq


      # Check on some odd edge cases in which checkpoints exist along with other behaviour.
      test_for_ckpt=$(/bin/ls -rt1 ../chk_point_dump/*${cnc}* 2> /dev/null | tail -1)
      if [ ! -z $test_for_ckpt ] ; then
        # Rare case where data has been completely written for the case, but the -unnecessary-
	# checkpoint is bad (that's the conclusion from an extant checkpoint and with an stdout
	# error: the error happened during the checkpoint write).
	# In this case, we'll consider it as a success for pragmatic reasons, since we have all 
	# of the data we require for the case.
	if grep -q 'Model run complete due to model time\|timestep completion' $lfile ; then
          echo -e "${YELLOW}  Badckpt! ~~ Kind of a ~~${GREEN} Success ${NC}"
            if [[ ! -z $act ]] ; then
              grep 'Miss \|un-enabled' $lfile > ${jn}.success
	      echo -e "\nWrite of final checkpoint, ${test_for_ckpt}, failed.\nSee ${lfile}\n" > ${jn}.success
            fi
          echo ""
          continue # check next case
        fi

        # Weird edge case where there is an error on cycle 1, but a checkpoint was created.
        #   In the strangest case, the checkpoint is good (this has happened).
        #   In the most normal case, the model failed during the checkpoint write.
        cycleid=$(sed 's/.*_//' <<< $lfile)
        if [ ${cycleid} == 1 ] ; then
          echo -e "${YELLOW}------------ WARNING ------------${NC}"
          echo "  There is an error on cycle 1, but a checkpoint was created."
          echo "  We will treat this as though the checkpoint is good and continue."
          echo "  However, it is likely that the checkpoint is bad, and the case"
          echo "  will fail on the start of the next cycle."
          echo -e "${YELLOW}------------ WARNING ------------${NC}"
          echo ""
          sleep 1
	  if [[ ! -z $act ]] ; then
            job_deletion ${jn}
            cd $owd
            if [ ${niq} -lt ${nqlim} ] ; then
              echo -e "${YELLOW}      ...resubmitting...${NC}"
              ${submit} ${inc}/tmp_qsub_dir/submonc_${cnc}
            else
              echo -e "${YELLOW}  -- Queue limit exceeded.  Job not submitted. --  ${NC}"
            fi
            cd $cwd
          fi # resubmit action
          continue  # check next case
        fi # cycle number 1 check
      fi # ckpt check

      # Look for previous restarts on this cycle
      aaa=0
      echo Searching ${jn}.rerun for "${lfile}".
      [ -f ${jn}.rerun ] && aaa=$(grep "${lfile}" ${jn}.rerun | wc -l)
      echo "  Found ${aaa} previous restarts"

      # If there are more than ${rtest} restart attempts, stop submitting
      if [[ $aaa -gt $rtest ]]; then
        echo " ...FAILED...that's enough"
        if [[ $system == "pbs" ]] ; then
          core=`grep 'ATP an' $lfile | cut -f 2 -d ' '`
          echo "gdb ${cwd}/../build_dir/monc_driver.exe core.atp.${core}."
        fi
echo -e "${RED}ATTENTION${NC}"
        echo ""
        continue # check next case
      fi

      # If there is already at least one previous restart, delete a checkpoint
      if [[ $aaa -ge 1 ]] && [[ $aaa -le $rtest ]] ; then
echo -e "${RED}ATTENTION${NC}"
        echo "    ... backing up one checkpoint (if present) ... "
        badckpt=$(/bin/ls -rt1 ../chk_point_dump/*${cnc}* 2> /dev/null | tail -1)
        echo "  Removing: $badckpt"
        [[ ! -z $act ]] && rm $badckpt
	# Remove entries in the rerun tracking file for the failed cycle so this 
	# process isn't accidentally repeated later.
	if [[ ! -z $act ]] ; then
	  grep -v "{lfile}" ${jn}.rerun > ${jn}.rerun.tmp
	  mv ${jn}.rerun.tmp ${jn}.rerun
	fi
        backup_flag=1
      fi

      # If one or none are found, save the output
      echo -e "${YELLOW}Recording ${lfile} in ${jn}.rerun."
      echo -e "      ...resubmitting...${NC}"
echo -e "${RED}ATTENTION. See:${NC}"
      if [[ ! -z $act ]] ; then
        echo ${lfile} >> ${jn}.rerun
        mv ${lfile} retain.${lfile}.${aaa}
        /bin/ls -alh retain.${lfile}.${aaa}

        # In the case where we step back a checkpoint, it's a good idea to also retain the previous stdout,
        #   as it is likely to also have failed, and it will be rerun anyay.
        if [[ ! -z $backup_flag ]] ; then
          lfilem1=$(/bin/ls -rt1 output*$cnc* 2> /dev/null | tail -1)
          [[ -f ${lfilem1} ]] && mv ${lfilem1} retain.${lfilem1}.badckpt
        fi

        # and resubmit
        job_deletion ${jn}
        cd $owd
        if [ ${niq} -lt ${nqlim} ] ; then
          ${submit} ${inc}/tmp_qsub_dir/submonc_${cnc}
        else
          echo -e "${YELLOW}  -- Queue limit exceeded.  Job not submitted. --  ${NC}"
        fi
        cd $cwd
      fi
      echo ""
      continue # check next case


    #----------------------------------------------------------------------------------  
    else # no errors were found

      if grep -qi 'Miss ' $lfile ; then
        echo ''
echo -e "${RED}ATTENTION${NC}"
        echo '  --- Miss match detected ---  '
echo -e "${RED}ATTENTION${NC}"
        echo ''
      fi

      if grep -q 'un-enabled' $lfile ; then
        echo ''
echo -e "${RED}ATTENTION${NC}"
        echo '  --- un-enabled detected ---  '
echo -e "${RED}ATTENTION${NC}"
        echo ''
      fi

      # Check for success (with no errors)  
      if grep -q 'Model run complete due to model time\|timestep completion' $lfile  ; then
      echo -e "${GREEN}  Success ${NC}"
        if [[ ! -z $act ]] ; then
          grep 'Miss \|un-enabled' $lfile > ${jn}.success
        fi
        echo ""
        continue # check next case
      fi

      # Check on the held job situation
      # Is a held job present?
      if eval $qstat | grep -w ${jn} | grep -q "${held_cond}" ; then
        echo "      ... only a held job remains ...restarting"
        echo " Progress checking stopped."
        echo " We do not yet understand this state of results."
        echo " However, it's possible that there is a job in the ERROR state,"
        echo " which is pretty ephemeral."
        echo " Take advantage of this situation, and investigate a solution."
        echo ''
exit 1
        if [[ ! -z $act ]] ; then
          job_deletion ${jn}
          cd $owd
          if [ ${niq} -lt ${nqlim} ] ; then
            ${submit} ${inc}/tmp_qsub_dir/submonc_${cnc}
          else
            echo -e "${YELLOW}  -- Queue limit exceeded.  Job not submitted. --  ${NC}"
          fi
          cd $cwd
        fi
        echo ""
        continue # check next case

      else # In this case, the held job failed, but we don't have an error message.
        echo "      ... held job failed...restarting"
        echo "See ../pbs_dir/$nextjob.OU for more info" 
        grep -B 1 '= PBS epilogue =' ../pbs_dir/$nextjob.OU | head -1
echo -e "${RED}ATTENTION${NC}"
        if [[ ! -z $act ]] ; then
          cd $owd
          if [ ${niq} -lt ${nqlim} ] ; then
            ${submit} ${inc}/tmp_qsub_dir/submonc_${cnc}
          else
            echo -e "${YELLOW}  -- Queue limit exceeded.  Job not submitted. --  ${NC}"
          fi
          cd $cwd
        fi
        echo ""
        continue # check next case
      fi # check for held jobs

    fi # end of check for errors and jobless states with out errors

  done # cnc loop over test cases

  cd $owd
  echo ''
  echo ''
  echo ''

done # inc loop over test harness configurations

exit 0
