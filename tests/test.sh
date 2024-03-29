#execute this from the root of the project `bash src/run.sh`
#(cd wdl-bamstats/&& ml cromwell/56-Java-11 && java -jar $EBROOTCROMWELL/cromwell.jar run  Dockstore.wdl -i test.wdl.json --workflow-root ./ )
#when running with slurm as a backend
#(cd wdl-bamstats/ && ml cromwell/56-Java-11 && java -Dconfig.file=../siteconfig/peregrine/slurm.conf -jar $EBROOTCROMWELL/cromwell.jar run Dockstore.wdl -i test.wdl.json --workflow-root ./ )
set -ex
(ml cromwell/56-Java-11 ;java -Dconfig.file=site/gearshift/slurm.conf -Xmx6g -jar $EBROOTCROMWELL/womtool.jar validate miR.wdl -i tests/integration/input.json)

(ml cromwell/56-Java-11; java -Dconfig.file=site/gearshift/slurm.conf -Xmx6g -jar $EBROOTCROMWELL/cromwell.jar run miR.wdl -i tests/integration/input.json --workflow-root ./)
touch tests/test.done
