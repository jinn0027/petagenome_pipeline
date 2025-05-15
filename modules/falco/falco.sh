#!/bin/bash

MY_DIR=$(cd $(dirname $BASH_SOURCE); pwd)

apptainer exec falco.sif falco -h





