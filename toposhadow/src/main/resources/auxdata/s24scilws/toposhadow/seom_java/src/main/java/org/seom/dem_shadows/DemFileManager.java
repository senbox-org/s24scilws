package org.seom.dem_shadows;

import org.orekit.rugged.errors.RuggedException;

public abstract class DemFileManager {

    /**
     * DEM root directory, should contains raster files
     */
    protected String demRootDir;

    /**
     * Default constructor
     * @param demRootDir dem root directory
     */
    public DemFileManager(String demRootDir) {
        this.demRootDir = demRootDir;
    }

    /**
     * check if there is at least one raster file in DEM root directory
     * @return true if at least one raster file is found in demRootDir
     * @throws S2GeolibException if no raster is found
     */
    public boolean findRasterFile() throws RuggedException {
        return findRasterFile(demRootDir);
    }

    /**
     * check if there is at least one raster file in the given directory
     * @param directory should contains raster files
     * @return true if at least one raster file is found
     * @throws S2GeolibException
     */
    protected abstract boolean findRasterFile(String directory) throws RuggedException;

    /**
     * Get a raster file path corresponding to given lat/lon values
     * @param latitude latitude (rad)
     * @param longitude longitude (rad)
     * @return the raster file path corresponding to given lat/lon values
     */
    protected abstract String getRasterFilePath(double latitude, double longitude);

    /**
     * To check if the DEM root dir is not an empty string
     * @return the DEM root directory
     */
    public String getDemRootDir(){
        return this.demRootDir;
    }
}
