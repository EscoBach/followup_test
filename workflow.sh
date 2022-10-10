#!/bin/bash

git fetch

BRANCH="$(git rev-parse --abbrev-ref HEAD)"

if test "$BRANCH" != "main"
then

	#create log file
	#log=logs/log_$(date +"%d-%m-%Y").txt
	#printf "Log File - " > $log
	#date >> $log

	R CMD BATCH main.r 

	git add .

	git commit -m 'auto_git'

	git push
fi


