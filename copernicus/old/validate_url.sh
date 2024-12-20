#!/bin/bash

# simple function to check http response code before downloading a remote file
# example usage:
# if `validate_url $url >/dev/null`; then dosomething; else echo "does not exist"; fi


function validate_url(){
  if [[ `wget --user ${2} --password ${3} -S --spider ${1} | grep 'tal'` ]]
  #if [[ `wget --user ${2} --password ${3} -S --spider ${1} | grep 'no such file'` ]]  
  then echo "true"
  fi
}
