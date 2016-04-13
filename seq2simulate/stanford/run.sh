#!/bin/bash

cd "$(dirname "$0")" && java -Xmx512m -cp sierra.jar:asi_fstrf.jar:commons-cli-1.1.jar:dom4j-1.6.1.jar:commons-collections-3.2.1.jar:org.apache.commons.lang.jar:jaxen-1.1.1.jar edu.stanford.hiv.utilities.createReport -k B0VQ-C341-5FPU-LO2P -f $1 -rf $2 -o "XMLVER1|QA ANALYSIS"
