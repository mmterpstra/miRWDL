#execute this from the root of the project `bash src/run.sh`
#(cd wdl-bamstats/&& ml cromwell/56-Java-11 && java -jar $EBROOTCROMWELL/cromwell.jar run  Dockstore.wdl -i test.wdl.json --workflow-root ./ )
#when running with slurm as a backend
#(cd wdl-bamstats/ && ml cromwell/56-Java-11 && java -Dconfig.file=../siteconfig/peregrine/slurm.conf -jar $EBROOTCROMWELL/cromwell.jar run Dockstore.wdl -i test.wdl.json --workflow-root ./ )

#smart defaults
PIPELINE="miR.wdl"
WORKFLOWROOT="$(realpath "$PWD")" 
CONFIG=$(ls $WORKFLOWROOT/site/hab*/cromwell.conf || ls $WORKFLOWROOT/site/*/cromwell.conf) 

while getopts i:s:w:r:f:d:p:c: flag
do
    case "${flag}" in
        #generic setting json
        i) INPUTJSON="$(realpath "${OPTARG}")" ;;
        #samplesheet json
        s) SAMPLEJSON="$(realpath "${OPTARG}")";;
        # usually the git dir so ./
        w) WORKFLOWROOT="$(realpath "${OPTARG}")";;
        # the location to store the (intermediate) files in and it should be accesable to the nodes it is being ran on (eg scratch)
        r) RUNROOT="$(realpath "${OPTARG}")";;
        # Location to be symlinked to the runroot for path abstraction
        f) FASTQRAWDIRS="$(realpath "${OPTARG}")";;
        # location of the data excluding the /data at the end also for path abstraction. So '/path/to/datadir/' should contain './data/'
        d) DATADIR="$(realpath "${OPTARG}")";;
        # omit to default to "miRWDL.wdl"
        p) PIPELINE=${OPTARG};;
        # Pipeline config override 
        c) CONFIG="$(realpath "${OPTARG}")";;
    esac
done

ml purge
(
    #init env
    
    set -e
    #create limited acces raw folder
    mkdir -p $RUNROOT/raw

    #init data
    if [ ! -e "$RUNROOT/data" ]; then
        #ln -s /groups/path/to/apps/data $RUNROOT/
        ln -s $DATADIR/data $RUNROOT/
        
    fi

    #archive inputs.json
    if [ ! -e "$RUNROOT/inputs.json" ]; then
        cp $INPUTJSON $RUNROOT/inputs.json
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON)
    elif [ -e "$RUNROOT/inputs.json" ]; then
        mv $RUNROOT/inputs.json{,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $INPUTJSON $RUNROOT/inputs.json
        echo "## "$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g')"## $PWD $0 $@" >> $SAMPLESHEETFOLDER/${HOSTNAME}_bestie.log
        
        #echo $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON){,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        mv $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON){,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $INPUTJSON)
    fi

    #archive sample.json
    if [ ! -e "$RUNROOT/sample.json" ]; then
        cp $SAMPLEJSON $RUNROOT/sample.json
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $SAMPLEJSON)
    elif [ -e "$RUNROOT/sample.json" ]; then
        mv $RUNROOT/sample.json{,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $SAMPLEJSON $RUNROOT/sample.json
        mv $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $SAMPLEJSON){,$(date --iso-8601=min | perl -wpe 's/[\+\:]/_/g').bak}
        cp $INPUTJSON $SAMPLESHEETFOLDER/${HOSTNAME}_$(basename $SAMPLEJSON)
    fi
    
    for FASTQRAWDIR in $(echo $FASTQRAWDIRS | tr , \ ); do
        if [ ! -e "$RUNROOT/raw/$(basename $FASTQRAWDIR)" ]; then
            ln -s $FASTQRAWDIR $RUNROOT/raw
        else 
            unlink $RUNROOT/raw/$(basename $FASTQRAWDIR)
            ln -s $FASTQRAWDIR $RUNROOT/raw
        fi
        # Now this should work for relative and complete paths for compatibility
        perl -i$(basename $FASTQRAWDIR).bak -wpe 's?'$FASTQRAWDIR'|./raw/'$(basename $FASTQRAWDIR)'?'$RUNROOT/raw/$(basename $FASTQRAWDIR)/'?g' $RUNROOT/sample.json
        #perl -i$(basename $FASTQRAWDIR).bak -wpe 's?'$FASTQRAWDIR'?'$RUNROOT/raw/$(basename $FASTQRAWDIR)/'?g' $RUNROOT/sample.json
    done
    #fix samplejson link?
    #ls -alh $RUNROOT
    perl -i.bak -wpe 's?.sampleJson": ".*",?.sampleJson": "'$RUNROOT/sample.json'",?g' $RUNROOT/inputs.json

    #start workflow
    ml cromwell/79-Java-11|| ml cromwell
    
    set -x
    java -Xmx8g -Dconfig.file=$CONFIG \
        -jar $EBROOTCROMWELL/womtool.jar validate \
        $WORKFLOWROOT/$PIPELINE \
        -i $RUNROOT/inputs.json 
    cd $RUNROOT
    
    echo "bash ${SCRIPTCALL}," $(git log -1 || echo "$(pwd);$(date)" ) | head -n 1  >> "$RUNROOT/nohup_"$(date --rfc-3339=date)".out"

    nohup java -Xmx8g -Dconfig.file=$CONFIG \
        -jar $EBROOTCROMWELL/cromwell.jar run \
        $WORKFLOWROOT/$PIPELINE \
        -i $RUNROOT/inputs.json \
        --workflow-root $WORKFLOWROOT &>> "$RUNROOT/nohup_"$(date --rfc-3339=date)".out"
    echo "watch log in $RUNROOT/nohup_"$(date --rfc-3339=date)".out"
)