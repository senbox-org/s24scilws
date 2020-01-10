package org.seom.dem_shadows;

import java.io.File;

//import org.apache.logging.log4j.LogManager;
//import org.apache.logging.log4j.Logger;
import org.gdal.gdal.Band;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconst;
import org.hipparchus.exception.DummyLocalizable;
import org.hipparchus.util.FastMath;
import org.orekit.rugged.errors.RuggedException;
import org.orekit.rugged.raster.UpdatableTile;


/**
 * Tile updater reading raster file using GDAL
 */
public class DemManager extends DEMTileUpdater {


    /**
     * File manager used to get raster file
     */
    private DemFileManager demFileManager;
    /**
     * Nb times the same lat/lon value is used to ask to update a tile (used to avoid infinite loop)
     */
    private int nbCall = 0;
    /**
     * Last lat value used to ask to update a tile
     */
    private double lastLat = Double.NaN;
    /**
     * Last lon value used to ask to update a tile
     */
    private double lastLon = Double.NaN;

    /**
     * Unique DemManager instance
     */
    private static DemManager singleton = null;

    /**
     * Get the unique instance
     * @return the unique instance of DemManager
     */
    public static synchronized DemManager getInstance() {
        return singleton;
    }

    /**
     * Init DemManager instance
     * @param demFileManager DEM file manager
     */
    public static synchronized DemManager initDemManager(DemFileManager demFileManager) {
        singleton = new DemManager(demFileManager);
        return singleton;
    }

    /**
     * Private constructor
     * @param demFileManager DEM file manager
     */
    private DemManager(DemFileManager demFileManager) {
        super();
        this.demFileManager = demFileManager;
    }



    /** Update given tile using DEM elevation
     * @param latitude latitude that must be covered by the tile (rad)
     * @param longitude longitude that must be covered by the tile (rad)
     * @param tile to update
     * @exception RuggedException if tile cannot be updated
     */
    @Override
    public void updateTile(double latitude, double longitude, UpdatableTile tile) throws RuggedException {

        // Check if latitude and longitude already called ...
        if (latitude == lastLat && longitude == lastLon) {
            nbCall++;
        } else {
            lastLat = latitude;
            lastLon = longitude;
            nbCall = 0;
        }
        if (nbCall > Constants.NB_TIME_ASK_SAME_TILE) {
            String str = String.format("infinite lopp for %3.8f long %3.8f lat ",longitude, latitude);
            DummyLocalizable message = new DummyLocalizable(str);
            throw new RuggedException(message);

        }

        String demFileName = "";
        try {
            demFileName = demFileManager.getRasterFilePath(latitude, longitude);
            File demFile = new File(demFileName);

            if (!demFile.exists()) {
                System.out.format("dem filename %s not found %n",demFileName);
                // No DEM => we are on water => read geoid
                GeoidManager geoidManager = GeoidManager.getInstance();
                geoidManager.updateTile(latitude, longitude, tile);
                return;
            }

            Dataset poDataset = gdal.Open(demFileName, gdalconst.GA_ReadOnly);

            int rasterLonSize = poDataset.GetRasterXSize();
            int rasterLatSize = poDataset.GetRasterYSize();
            double[] geoTransformInfo = poDataset.GetGeoTransform();

            double minLat = Double.NaN;
            double minLon = Double.NaN;
            double lonStep = Double.NaN;
            double latStep = Double.NaN;
            if (geoTransformInfo.length < 6) { // Check that the geoTransformInfo is correct
                String str = "GDALGeoTransform does not contains 6 elements";
                DummyLocalizable message = new DummyLocalizable(str);
                throw new RuggedException(message);

            } else {
                // Get the step and range in latitude and longitude
                lonStep = FastMath.abs(geoTransformInfo[1]);
                latStep = FastMath.abs(geoTransformInfo[5]);

                minLon = geoTransformInfo[0] + 0.5 * lonStep;
                minLat = geoTransformInfo[3] - latStep * (rasterLatSize - 0.5);
            }

            Band band = poDataset.GetRasterBand(1);

            // Init Tile Geometry
            tile.setGeometry(FastMath.toRadians(minLat), FastMath.toRadians(minLon), FastMath.toRadians(latStep), FastMath.toRadians(lonStep), rasterLatSize, rasterLonSize);

            // Loop over dem values
            double[] data = new double[rasterLatSize * rasterLonSize];

            // Check if the DEM is over the geoid or not
            GeoTransformTools gtTools = null;
            GeoidManager geoidManager = null;
            boolean isDemRelativeToGeoid = isDemRelativeToGeoid();
            if (isDemRelativeToGeoid) { // SRTM is given above the geoid
                geoidManager = GeoidManager.getInstance();
                gtTools = new GeoTransformTools(geoTransformInfo, rasterLonSize, rasterLatSize);
            }

            // We read all raster at once : we try to limit the numbers of Band.ReadRaster calls
            // because when too much jvms are launch together and try to read some raster,
            // performances are not good at all
            band.ReadRaster(0, 0, rasterLonSize, rasterLatSize, data);

            // Get the no data value from the raster
            Double[] noDataValue = new Double[1];
            band.GetNoDataValue(noDataValue);

            // test if the no data value exists
            Boolean noDataValueExist = false;
            if (noDataValue[0] != null ) noDataValueExist = true;

            // from bottom left to upper right corner
            for (int iLat = rasterLatSize - 1; iLat >= 0; iLat--) {

                for (int jLon = 0; jLon < rasterLonSize; jLon++) {

                    double elevationOverEllipsoid = 0.0;
                    elevationOverEllipsoid = data[iLat * rasterLonSize + jLon];

                    if (noDataValueExist && (elevationOverEllipsoid == noDataValue[0])){
                        elevationOverEllipsoid = 0.0;
                    }
                    if (isDemRelativeToGeoid) {
                        // The elevation value we send to rugged must be computed against ellipsoid
                        // => when DEM is SRTM , we must add geoid value
                        double lon = gtTools.getXFromPixelLine(jLon, iLat);
                        double lat = gtTools.getYFromPixelLine(jLon, iLat);
                        elevationOverEllipsoid = elevationOverEllipsoid + geoidManager.getElevation(lat, lon);
                    }
                    // Set elevation over the ellipsoid
                    tile.setElevation(rasterLatSize - 1 - iLat, jLon, elevationOverEllipsoid);
                }
            }
            band.delete();
            band = null;
            poDataset.delete();
            poDataset = null;

        } catch (RuggedException e) {
            throw e;
        }
    }

    /**
     * If the DEM is given relative to the geoid (for instance SRTM)
     * @return true if DEM is relative to Geoid (Case of SRTM); false in other cases
     */
    public boolean isDemRelativeToGeoid() {
        boolean returned = false;
        if (demFileManager instanceof SRTM1sFileManager || demFileManager instanceof RasterFileManager) {
            returned = true;
        }
        return returned;
    }



}
