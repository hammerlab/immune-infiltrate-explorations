#!/bin/bash

export CMDSTAN_PATH=/home/elizachang/bin/cmdstan-2.14.0

export MODEL_DIR=/home/elizachang/immune-infiltrate-explorations/model-single-origin-samples/models
export MODEL_NAME=$1
export MODEL_EXE_PATH=$MODEL_DIR/$MODEL_NAME
export DATA_PATH=$2

cd $CMDSTAN_PATH && make $MODEL_EXE_PATH

cd $MODEL_DIR

./$MODEL_NAME sample data file=$DATA_PATH
