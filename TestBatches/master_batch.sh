#!/bin/bash

IGNORE_TESTS=()
IGNORE_TESTS+=(chimera.bellerophon)

TEST_DIR=TestBatches
if [ "x$1" != "x" ] ; then
	TEST_DIR=$1
fi

MOTHUR_EXEC=mothur

for TEST_FILE in `find $TEST_DIR -type f -name batch` ; do
        PROCESS=True
	for IGNORE_TEST in ${IGNORE_TESTS[@]}; do
		if [ `echo $TEST_FILE | grep -v /$IGNORE_TEST/batch | wc -l` == 0 ] ; then
			PROCESS=False
		fi
	done
	if [ ${PROCESS} == 'True' ] ; then
		echo "Processing $TEST_FILE"
		./mothur "$TEST_FILE"
	else
		echo "Ignoring $TEST_FILE"
	fi
done
