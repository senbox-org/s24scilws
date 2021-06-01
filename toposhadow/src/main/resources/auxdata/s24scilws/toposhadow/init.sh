
# Path to OrthotoolPython
export PYTHONPATH=/usr/bin/python3.6
export SEOM_TOPO_DIR=${HOME}/dev/s24scilws/toposhadow/src/main/resources/auxdata/s24scilws/toposhadow

# Data
export SEOM_TOPO_DATA=${SEOM_TOPO_DIR}/seom_data/

# Add orthotool python script to PYTHONPATH
export SEOM_TOPO_PYTHON="${SEOM_TOPO_DIR}/seom_python/src/"
export PYTHONPATH="${SEOM_TOPO_DIR}:$PYTHONPATH"

export SEOM_JAVA="${SEOM_TOPO_DIR}/seom_java/target/seom-1.0-jar-with-dependencies.jar"
export GDAL_DIR=$HOME/custom_gdal/COTS

# GDAL binding's java library path
export GDAL_LIBRARY_PATH=$GDAL_DIR/gdal-2.4.0/swig/java