#!/bin/bash
# -*- coding:utf-8 -*-

for fname in $(ls $1/*.txt)
do
   starcode -d2 -t8 $fname > ${fname%txt}stc
done
