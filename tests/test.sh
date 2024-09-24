#execute this from the root of the project `bash src/run.sh`
#(cd wdl-bamstats/&& ml cromwell/56-Java-11 && java -jar $EBROOTCROMWELL/cromwell.jar run  Dockstore.wdl -i test.wdl.json --workflow-root ./ )
#when running with slurm as a backend
#(cd wdl-bamstats/ && ml cromwell/56-Java-11 && java -Dconfig.file=../siteconfig/peregrine/slurm.conf -jar $EBROOTCROMWELL/cromwell.jar run Dockstore.wdl -i test.wdl.json --workflow-root ./ )
set -ex
if [ "x$TMP01" != "x" ] ; then
    DATADIR=$TMP01/apps/
    TESTPREFIX="_gearshift"
elif [ "x$TMPSCRATCH" != "x" ]; then
    DATADIR=$TMPSCRATCH/apps/
    #ESTPREFIX="_habrok"
fi
bash run_project.sh -i tests/integration/input${TESTPREFIX}.json -s ./tests/integration/config/samples.json -r ./tests/runs/local/ -f ./tests/data/fastq -d $DATADIR
touch tests/test.done
