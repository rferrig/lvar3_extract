#!/bin/bash
echo "Grabbing sequences from $1"
bash extract_transcripts.sh $1;
bash get_100bp_up.sh $1;
bash get_100bp_down.sh $1;
python combine_100bp_up_and_transcript.py
