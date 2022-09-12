#!/usr/bin/env bash
HERE=$(cd "$(dirname "$0")" && pwd -P)
FILE=cp2k-perf.txt

rm -f $FILE
for MVALUE in 9 ; do
    for NVALUE in 9 ; do
        for KVALUE in 9 ; do
        echo "M = ${MVALUE}, N = ${NVALUE}, K = ${KVALUE}"
        "${HERE}/cp2k-tester" "${MVALUE}" "${SIZE}" 0 "${NVALUE}" "${KVALUE}" >>"${FILE}";  2>&1
        echo >>"${FILE}"
        done
    done  
done

