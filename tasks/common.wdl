version 1.0

task CollapseFastq {
    input {
        Array[File] reads
        #String dockerImage = "quay.io/biocontainers/cutadapt:4.2--py310h1425a21_0"
        String outputPrefix
        String? threeLetterName 
        Int memoryGb = 1
        Int timeMinutes = 10 + ceil(size(reads, "GiB") * 40 )
        String? dockerImage
        #String dockerImage = "quay.io/biocontainers/coreutils" #idk needs perl 5 and coreutils
    }
    command {
        set -e -o pipefail
        mkdir -p "$(dirname ${outputPrefix})"
        zcat ${sep=" " reads} | \
        perl -wne 'chomp; print $_."\n" if($.%4==2)' | \
        sort --temporary-directory="$(dirname ${outputPrefix})" | \
        uniq -c | \
        python -c "import sys;[sys.stdout.write('>'+sys.argv[1]+'_'+str(count)+'_x'+'\n'.join(line.split(' ')[-2:])) for count,line in enumerate(sys.stdin)]" AAA > \
         "${outputPrefix}"".md.fa"
    }
    output {
        File outputCollapsedFasta = outputPrefix + ".md.fa"
    }
        runtime {
        memory: memoryGb*1024
        timeMinutes: timeMinutes
        #ocker: dockerImage
    }

}


#Borrowed stuff below here

#https://github.com/biowdl/tasks/blob/c92755e510723da731ba92637c41e58c8490b5fc/common.wdl#L66
# Copyright (c) 2017 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

task AppendToStringArray {
    input {
        Array[String] array
        String string
        
        Int memory = 256
    }

    command {
        echo "~{sep='\n' array}
        ~{string}"
    }

    output {
        Array[String] outArray = read_lines(stdout())
    }

    runtime {
        memory: memory
    }
}

task CreateLink {
    # Making this of type File will create a link to the copy of the file in
    # the execution folder, instead of the actual file.
    # This cannot be propperly call-cached or used within a container.
    input {
        String inputFile
        String outputPath

        Int memory = 256
    }

    command {
        echo $PWD
        ln -sf "$(realpath ~{inputFile})"  "~{outputPath}"
    }

    output {
        File link = outputPath
    }

    runtime {
        memory: memory
    }
}

task ConcatenateTextFiles {
    input {
        Array[File] fileList
        String combinedFilePath
        Int memory = 256
    }

    # When input and output is both compressed decompression is not needed.
    #String cmdPrefix = if (unzip && !zip) then "zcat " else "cat "
    #String cmdSuffix = if (!unzip && zip) then " | gzip -c " else ""

    command {
        set -e -o pipefail
        mkdir -p "$(dirname ~{combinedFilePath})"
        cat ~{sep=" " fileList} > ~{combinedFilePath}
    }

    output {
        File combinedFile = combinedFilePath
    }

    runtime {
        memory: memory
    }
}