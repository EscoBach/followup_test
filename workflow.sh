#!/bin/bash

cd followup_test/

git fetch

COMMIT=$(git rev-list --right-only --count HEAD...origin/main)

if (($COMMIT > 0))
then

	#create log file
	#log=logs/log_$(date +"%d-%m-%Y").txt
	#printf "Log File - " > $log
	#date >> $log

	git pull

	R CMD BATCH main.r 

	git add .

	git commit -m 'auto_git'

	git push
fi


