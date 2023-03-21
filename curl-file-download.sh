#!/usr/bin/env bash
# curl-file-download.sh - Script to download files using curl

# AUTHOR
#
# Richard White <rich@buschlab.org>
#
# COPYRIGHT AND LICENSE
#
# This software is Copyright (c) 2023. Queen Mary University of London.
#
# This is free software, licensed under:
#
#  The GNU General Public License, Version 3, June 2007

# define function for errors
# arg1 exit code
# arg2 Message for success
# arg3 Message for failure
function error_checking (){
    if [ $1 -eq 0 ]
    then
        if [[ -n $2 && $verbose -eq 1 ]]; then echo $2 ; fi
    else
        date >& 2
        echo $3 >&2
        exit $1
    fi
}

USAGE="curl-file-download.sh [options] input_file

The input file should be a tab-separated file of urls and optionally file names"

OPTIONS="Options:
    -d    print debugging info
    -v    verbose output
    -h    print help info"

# default values
debug=0
verbose=0

while getopts ":dhv" opt; do
  case $opt in
    d)
      debug=1
      ;;
    h)
      echo ""
      echo "$USAGE"
      echo "$OPTIONS"
      exit 1
      ;;
    v)
      verbose=1
      ;;
    \?)
      echo ""
      echo "Invalid option: -$OPTARG" >&2
      echo "$USAGE" >&2
      echo "$OPTIONS" >&2
      exit 1
      ;;
    :)
      echo ""
      echo "Option -$OPTARG requires an argument!" >&2
      echo "$USAGE" >&2
      echo "$OPTIONS" >&2
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

# unpack arguments
INPUT_FILE=$1

while read line
do
    url=$( echo $line | awk '{ print $1 }' )
    filename=$( echo $line | awk '{ print $2 }' )
    if [[ -z $filename ]]; then
        filename=$( echo $url | sed -e 's|^.*/||')
    fi
    if [[ $debug -eq 1 ]]; then
        echo URL:$url FILE:$filename
        echo "CMD: curl -sS -o $filename.tmp $url && mv $filename.tmp $filename"
    fi
    curl -sS -o $filename.tmp $url && mv $filename.tmp $filename
    error_checking $? "Download of file, $filename, succeeded." "Download of file, $filename, failed: $?"
done < $INPUT_FILE
