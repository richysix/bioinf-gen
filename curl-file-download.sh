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

source bash_functions.sh

USAGE="curl-file-download.sh [options] input_file

The input file should be a tab-separated file of urls and optionally file names"

OPTIONS="Options:
    -d    print debugging info
    -v    verbose output
    -h    print help info"

# default values
debug=0
verbose=0

while getopts ":dhvq" opt; do
  case $opt in
    d) debug=1 ;;
    h)  usage "" ;;
    v)  verbose=1 ;;
    q)  verbose=0 ;;
    \?) usage "Invalid option: -$OPTARG" ;;
    :)  usage "Option -$OPTARG requires an argument!" ;;
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
    if [[ -e $filename ]]; then continue; fi
    CMD="curl -sS -o $filename.tmp $url && mv $filename.tmp $filename"
    if [[ $debug -eq 1 ]]; then
        echo URL:$url FILE:$filename
        echo "$CMD"
    fi
    eval $CMD
    SUCCESS=$?
    echo $SUCCESS
    SUCCESS_MSG="Download of file, $filename, succeeded."
    ERROR_MSG="Download of file, $filename, failed: $?"
    error_checking $SUCCESS
    sleep 2
done < $INPUT_FILE
