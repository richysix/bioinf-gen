#!/usr/bin/env bash
# bash_functions.sh - Script containing some common functions

# function to display usage and optional message
# arg1 usage statemnt
# arg2 options
# arg3 Message
function usage (){
    echo $3 >&2
    echo $1 >&2
    echo $2 >&2
    exit 2
}

# function to check the exit status of a command for errors
# arg1 exit code
# arg2 Message for success
# arg3 Message for failure
function error_checking (){
    if [[ $1 -eq 0 ]]
    then
        if [[ -n $2 && $verbose -ge 1 ]]; then echo $2 ; fi
    else
        date >& 2
        echo $3 >&2
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

