#!/bin/bash 

sort -k8,8g -k6,6g $1  | awk '{if ($8<=-2.0 || $8 == "-inf" ) print}' 
