#!/usr/bin/env bash
# bash_functions.sh - Script containing some common functions

# function to display usage and optional message
# arg1 Message
# Assumes variable USAGE and OPTIONS have been set
function usage (){
    echo $1 >&2
    echo -e "$USAGE\n" >&2
    echo "$OPTIONS" >&2
    exit 2
}

# function to check the exit status of a command for errors
# arg1 exit code
# arg2 Message for success
# arg3 Message for failure
function error_checking (){
    if [[ -z $SUCCESS_MSG ]]; then SUCCESS_MSG=$2; fi
    if [[ -z $ERROR_MSG ]]; then ERROR_MSG=$3; fi
    if [[ $1 -eq 0 ]]
    then
        if [[ -n $SUCCESS_MSG && $verbose -ge 1 ]]; then echo $SUCCESS_MSG ; fi
    else
        date >& 2
        echo $ERROR_MSG >&2
        exit $1
    fi
}

function file_checking (){
    if [ "$1" = "f" ]
    then
        if [[ ! -e $2 || ! -r $2 ]]; then
            error_checking 2 "" "$3: $2 either does not exist or is not readable!"
        fi
    else
        if [[ ! -d $2 || ! -x $2 ]]; then
            error_checking 2 "" "$3: $2 either does not exist or is not executable!"
        fi
    fi
}

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

