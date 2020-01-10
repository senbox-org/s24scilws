# SEOM Topography shadow

We propose a simple approach for detection and masking of topographic shadows. The principle of shadow identification is to map a ground area and estimate for each sample ground point whether it is directly under Sun light or if it is behind another topographic feature as seen from the Sun.

### Prerequisites

The programming language version required are:
- Python 2.7.5
- Java 1.8.0_101
- maven

Libraries and tools:
- Rugged 2.0 
- Orekit 9.1 
- Hipparchus 1.3
- GDAL : 2.4.0 (bindings java and python)

### Installing and downloading libraries

# GDAL

Get the binaries here:
https://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries

Binaries may not be available, you will have to compile the sources. You can follow the building process here: https://github.com/wri/forma-clj/wiki/Building-GDAL.

create the GDAL_DIR environnement variable (the location where the package will be installed)

```
cd $GDAL_DIR
wget http://download.osgeo.org/gdal/2.4.0/gdal-2.4.0.tar.gz
tar -xvzf gdal-2.4.0.tar.gz
cd gdal-2.4.0
./configure --with-java=yes
./make
```

then edit gdal-2.4.0/swig/java/java.opt with your JAVA_HOME path

```
cd gdal-2.4.0/swig/java/
./make
```
following environement variable  will be necessary for launching topograpgy shadow extraction

``
GDAL_LIBRARY_PATH = ${GDAL_DIR}/gdal-2.4.0/swig/java
``

# SEOM_JAVA

```
cd seom_java
mvn clean install
```

file pom.xml contains all dependancy informations

following environement variable  will be necessary for launching topograpgy shadow extraction

``
 SEOM_JAVA="seom_java/target/seom-1.0-jar-with-dependencies.jar"
``
note that jar are available in seom_java_jar direcotry

## Configuration

following data is necessary :
* orekit data : https://gitlab.orekit.org/orekit/orekit-data/-/archive/master/orekit-data-master.zip
* geoid path (EGM96) : 
* SRMT1s tiles : available here (https://earthexplorer.usgs.gov) 

following DEM structure shoould be respected for DEM data :
```
├── w001
│   ├── w001n49.tif
.
.
├── e001
│   ├── e001n49.tif
│   ├── e001n50.tif
│   ├── e001n51.tif
│   ├── e001n52.tif
.
.
├── e142
│   └── e142s39.tif
├── e143
│   └── e143s39.tif


```

if a tile is missing in processing, the a warning is emited


these elements have to be filled in :

```
seom_data/conf/rugged.conf
```

## Running 

You have to set the environment variables before launching the python script by sourcing the init.sh file. 
/!\ You will probably have to change the GDAL_DIR library path according to your gdal directory path.

This is the command to launch the seom topography mask extraction:
```
python ./seom_python/src/create_topo_mask.py

usage: create_topo_mask.py [-h] [-v] [--ULLR ULLR ULLR ULLR ULLR]
                           [--EPSG EPSG] [--threshold THRESHOLD]
                           ruggedConf product output tmp
                           
  ruggedConf            rugged configuration file
  product               product (.SAFE root dir) or L1C/L2A granule MTD (.xml)
  output                output mask filename (vector data in ESRI shaprefile
                        format)
  tmp                   working dir

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Verbose mode
  --ULLR ULLR ULLR ULLR ULLR
                        user defined footprint [ULX ULY LRX LRY]
  --EPSG EPSG           EPSG code of --ULLR
  --threshold THRESHOLD
                        detection threshold (0.1 meter by default)

```

## Example 

an example script is available in :

```
./seom_python/scripts/launch_seom_topo.sh
```


