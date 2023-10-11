# ProBound
---
This repository contains the source code for ProBound.

To compile ProBound, first install Maven and then execute:

`cd ProBoundApp`
`mvn package`

This will create the compiled JAVA file `ProBoundApp/target/ProBound-jar-with-dependencies.jar`

To run ProBound, first set the environmental variable

`export PROBOUND_DIR="/path/to/GitHub/ProBound/"`

Next, the configuration builder file `config.json` needs to be converted into a format that can be used by ProBound:

`java -jar $PROBOUND_DIR/ProBoundApp/target/ProBound-jar-with-dependencies.jar -b -c config.json > config.full.json`

Examples of the configuration builder and the converted file can be found under `sample_data`.

Finally, run ProBound using:

`java -jar $PROBOUND_DIR/ProBoundApp/target/ProBound-jar-with-dependencies.jar -c config.full.json`

