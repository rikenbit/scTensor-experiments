#!/bin/sh

ini=$1
fin=$2
for((i = ini; $i < $fin +1; i ++))
do
scancel $i
done 