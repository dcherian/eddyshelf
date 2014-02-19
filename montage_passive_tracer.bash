#!/bin/bash

DIR=videos/runew-03-close-z=
mkdir temp
montage $DIR"0"/$FNAME $DIR"-100"/$FNAME $DIR"-250"/$FNAME $DIR"-1000"/$FNAME -geometry 1600x900 output.jpg
