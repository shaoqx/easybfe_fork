#!/usr/bin/env bash

WDIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

run_step_seq() {
  local step_dir="$1"
  local name="$2"

  cd "$step_dir" || return 1

  echo "Running $step_dir ..."
  source "$name.sh" > "$name.stdout" 2>&1
  local rc=$?

  if [ $rc -ne 0 ]; then
    mv "$WDIR/running.tag" "$WDIR/error.tag"
    echo "Error occurs in $name (exit code $rc)"
    cd "$WDIR"
    return $rc
  fi

  cd "$WDIR"
}


cleanup() {
    echo "[`date`] Caught termination signal. Cleaning up..."
    if [ -f $WDIR/running.tag ]; then mv $WDIR/running.tag $WDIR/killed.tag; fi
    echo "[`date`] Cleanup done."
    exit 2 
}
trap cleanup TERM INT HUP

start=$(date +%s)
if [ -f running.tag ]; then
  echo "Found running.tag: a run may still be in progress. Skip this run."
  echo "Delete running.tag to force rerun."
  exit 0
fi
if [ -f done.tag ]; then
  echo "Found done.tag: previous run has finished. Skip this run."
  echo "Delete done.tag to force rerun."
  exit 0
fi
if [ -f error.tag ]; then
  echo "Found error.tag: previous run ended with error. Skip this run."
  echo "Delete error.tag to force rerun."
  exit 0
fi
touch running.tag

echo "Running 01.em"
run_step_seq lambda0/01.em 01.em || exit 1
run_step_seq lambda1/01.em 01.em || exit 1
run_step_seq lambda2/01.em 01.em || exit 1
run_step_seq lambda3/01.em 01.em || exit 1
run_step_seq lambda4/01.em 01.em || exit 1
run_step_seq lambda5/01.em 01.em || exit 1
run_step_seq lambda6/01.em 01.em || exit 1
run_step_seq lambda7/01.em 01.em || exit 1
run_step_seq lambda8/01.em 01.em || exit 1
run_step_seq lambda9/01.em 01.em || exit 1
run_step_seq lambda10/01.em 01.em || exit 1
run_step_seq lambda11/01.em 01.em || exit 1
run_step_seq lambda12/01.em 01.em || exit 1
run_step_seq lambda13/01.em 01.em || exit 1
run_step_seq lambda14/01.em 01.em || exit 1
run_step_seq lambda15/01.em 01.em || exit 1

echo "Running 02.heat"
mpirun -np 16 pmemd.cuda.MPI -ng 16 -groupfile 02.heat.groupfile
if [ $? -ne 0 ]; then
  mv running.tag error.tag && echo "Error occurs!"
  exit 1
fi


echo "Running 03.pres"
mpirun -np 16 pmemd.cuda.MPI -ng 16 -groupfile 03.pres.groupfile
if [ $? -ne 0 ]; then
  mv running.tag error.tag && echo "Error occurs!"
  exit 1
fi


echo "Running 04.pre_prod"
mpirun -np 16 pmemd.cuda.MPI -ng 16 -groupfile 04.pre_prod.groupfile
if [ $? -ne 0 ]; then
  mv running.tag error.tag && echo "Error occurs!"
  exit 1
fi


echo "Running 05.prod"
mpirun -np 16 pmemd.cuda.MPI -ng 16 -groupfile 05.prod.groupfile -rem 3 -remlog 05.prod.log
if [ $? -ne 0 ]; then
  mv running.tag error.tag && echo "Error occurs!"
  exit 1
fi



mv running.tag done.tag

end=$(date +%s)
duration=$((end - start))

hours=$(( duration / 3600 ))
minutes=$(( (duration % 3600) / 60 ))
seconds=$(( duration % 60 ))

echo "Execution time: ${hours} h ${minutes} min ${seconds} sec"