#!/bin/sh

###################################################################
# Wrapper to run the regression tests for the assembler
###################################################################
CC=$(which cc)

BINDIR=""
DIRBIN_LIST="AS_GKP:gatekeeper AS_CNS:consensus"

is_64bit()
{
  local ret;
  local src_file=/tmp/64$$.c
  local exe_file=/tmp/64$$.exec

  cat <<EOF > $src_file
int main() { return (sizeof(unsigned long) * 8) != 64; }
EOF

  $CC -o $exe_file $src_file > /dev/null 2>&1
  if [ ! -f $exe_file ]; then
    echo "is_64bit: couldn't compile test program"
    exit 1
  fi

  $exe_file
  ret=$?

  rm -f $src_file
  rm -f $exe_file

  return $ret
}

is_32bit()
{
  local ret;
  local src_file=/tmp/32$$.c
  local exe_file=/tmp/32$$.exec

  cat <<EOF > $src_file
int main() { return (sizeof(unsigned long) * 8) != 32; }
EOF

  $CC -o $exe_file $src_file > /dev/null 2>&1
  if [ ! -f $exe_file ]; then
    echo "is_32bit: couldn't compile test program"
    exit 1
  fi

  $exe_file
  ret=$?

  rm -f $src_file
  rm -f $exe_file

  return $ret
}

set_BINDIR() {
  if is_64bit ; then 
    pushd ../../Linux64/bin > /dev/null 2>&1
    BINDIR=$(pwd)
    popd > /dev/null 2>&1
  elif is_32bit ; then
    pushd ../../Linux/bin   > /dev/null 2>&1
    BINDIR=$(pwd)
    popd > /dev/null 2>&1
  else 
    echo "Hrm: can't tell if this is 32 or 64 bit"
    exit 1
  fi
}

check_BINS() {
  local dirbin=""
  local bin=""

  for dirbin in $DIRBIN_LIST
  do
    bin=${BINDIR}/$(echo $dirbin | awk -F: '{ print $2 }') 

    if [ ! -f $bin ]; then
      echo "check_BINS: file $bin does not exist"
      exit 1
    fi

    if [ ! -x $bin ]; then
      echo "check_BINS: file $bin is not executable"
      exit 1
    fi
  done
}

####
# MAIN:
####
set_BINDIR
check_BINS

for dirbin in $DIRBIN_LIST
do
  dir=$(echo $dirbin | awk -F: '{ print $1 }') 
  bin=${BINDIR}/$(echo $dirbin | awk -F: '{ print $2 }') 

  if [ -d $dir ]; then
    pushd $dir > /dev/null 2>&1
    for test in $(ls -1 [0-9][0-9][0-9].sh)
    do
      echo "[*] Running $dir/$test ..."
      cwd=$(pwd)
      ${cwd}/$test $bin
      test_status=$?
      if [ $test_status != "0" ]; then
        echo "Test $test failed"
      fi
    done
    popd > /dev/null 2>&1
  else
    echo "Error: directory $dir does not exist"
  fi
done

