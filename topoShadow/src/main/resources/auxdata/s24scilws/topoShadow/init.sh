
# Path to OrthotoolPython
export SEOM_TOPO_DIR=$PWD

# Data
export SEOM_TOPO_DATA=${SEOM_TOPO_DIR}/seom_data/


# Add orthotool python script to PYTHONPATH
export SEOM_TOPO_PYTHON="${ORTHOTOOL_DIR}/seom_python/src/"
export PYTHONPATH="${SEOM_TOPO_DIR}:$PYTHONPATH"

export SEOM_JAVA="${SEOM_TOPO_DIR}/seom_java/target/seom-1.0-jar-with-dependencies.jar"

# GDAL binding's java library path
export GDAL_LIBRARY_PATH="${GDAL_DIR}/COTS/gdal-2.4.0/swig/java"

