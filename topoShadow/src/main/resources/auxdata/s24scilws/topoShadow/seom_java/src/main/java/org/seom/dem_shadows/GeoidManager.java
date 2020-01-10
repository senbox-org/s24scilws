package org.seom.dem_shadows;

import org.gdal.gdal.Band;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconst;
import org.hipparchus.exception.DummyLocalizable;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.rugged.errors.RuggedException;
import org.orekit.rugged.raster.UpdatableTile;



/**
 * Tile updater reading the Geoid using GDAL.
 * <pre>
 * Availability of the tiles:
 *     [ -180 deg, +180 deg [ in longitude
 *     [  -90 deg, + 90 deg [ in latitude
 * (See Tile.Location.HAS_INTERPOLATION_NEIGHBORS constraints)
 * </pre>
 * Use the lazy evaluation design pattern.
 *
 * @author Guylaine Prat
 */
public class GeoidManager extends DEMTileUpdater {

    /**
     * Geoid file path
     */
    private String geoidFilePath;

    /**
     * Geoid elevations for the whole Earth (read from Geoid dataset)
     * Filled in with GDAL convention : [0][0] = North West
     */
    private double[][] elevationsWholeEarth = null;

    /**
     * Step (> 0) in longitude (deg)
     */
    private double lonStepInDeg = Double.NaN;

    /**
     * Step (> 0) in latitude (deg)
     */
    private double latStepInDeg = Double.NaN;

    /**
     * Geoid sub-tiles definition
     */
    private GeoidTileDefinition tile_NORTH_FARWEST = null;
    private GeoidTileDefinition tile_NORTH_WEST = null;
    private GeoidTileDefinition tile_NORTH_EAST = null;
    private GeoidTileDefinition tile_NORTH_FAREAST = null;
    private GeoidTileDefinition tile_SOUTH_FARWEST = null;
    private GeoidTileDefinition tile_SOUTH_WEST = null;
    private GeoidTileDefinition tile_SOUTH_EAST = null;
    private GeoidTileDefinition tile_SOUTH_FAREAST = null;

    /**
     * Earth quadrants limits in latitude and longitude
     */
    private EarthQuadrantsLimits earthQuadrantsLimits = null;

    /**
     * Size in longitude of the Geoid (number of columns/pixels)
     */
    private int rasterLonSize;

    /**
     * Size in latitude of the Geoid (number of lines)
     */
    private int rasterLatSize;

    /**
     * Unique GeoidManager instance
     */
    private static GeoidManager singleton = null;

    /**
     * Get the unique instance
     * @return the unique instance of GeoidManager
     */
    public static synchronized GeoidManager getInstance() {
        return singleton;
    }

    /**
     * Init GeoidManager instance
     * @param geoidFilePath geoid file path
     */
    public static synchronized GeoidManager initGeoidManager(String geoidFilePath) {
        singleton = new GeoidManager(geoidFilePath);
        return singleton;
    }

    /**
     * Private constructor
     * @param geoidFilePath geoid file path
     */
    private GeoidManager(String geoidFilePath) {
        super();
        this.geoidFilePath = geoidFilePath;
    }

    /**
     * Update a given tile using Geoid elevation
     * @param latitude (rad)
     * @param longitude (rad)
     * @param tile tile to update
     * @throws RuggedException if tile cannot be updated
     */
    @Override
    public void updateTile(double latitude, double longitude, UpdatableTile tile) throws RuggedException {

        try {
            // Normalize latitude between ±π/2 and longitude between ±π
            // (needed for the algorithms below)
            double latitudeNorm = FastMath.toDegrees(MathUtils.normalizeAngle(latitude, FastMath.PI / 2));
            double longitudeNorm = FastMath.toDegrees(MathUtils.normalizeAngle(longitude, 0));

            // Read the whole Geoid and split the Earth in 8 quadrants (for Rugged tiles cache)
            if (this.elevationsWholeEarth == null) { // done only if necessary at the first call (lazy evaluation design pattern)
                initGeoidSubTiles();
            } // end if this.elevationsWholeEarth == null

            // Search in which quadrant the (latitude, longitude) belongs to
            EarthQuadrant currentQuadrant = getQuadrant(latitudeNorm, longitudeNorm);

            // Identify the Geoid sub-tile the (latitude, longitude) belongs to
            GeoidTileDefinition foundGeoidTile = identifyGeoidTile(currentQuadrant);

            // Fill in the tile
            // ================
            // Fill the tile and get the max index in longitude when filling the tile: needed for FarEast tiles ...
            int indexLonMax = fillTheTile(foundGeoidTile, tile);

            // Special cases for the FarEast Tile
            // ----------------------------------
            // For FarEast tiles, one extra column was added to the longitude number of column
            // that must be filled in with the first column of the FarWest corresponding tile
            Boolean isFarEast = currentQuadrant.isFarEast();
            if (isFarEast){
                EarthHemisphere currentHemisphere = currentQuadrant.getHemisphereNS();
                if (currentHemisphere == EarthHemisphere.NORTH){
                    GeoidTileDefinition northFarWestTile = this.tile_NORTH_FARWEST;
                    fillTheLastColumnTile(northFarWestTile, tile, indexLonMax);
                } else { // currentHemisphere == EarthHemisphere.SOUTH
                    GeoidTileDefinition southFarWestTile = this.tile_SOUTH_FARWEST;
                    fillTheLastColumnTile(southFarWestTile, tile, indexLonMax);
                }
            }

        } catch (RuggedException e) { // Must convert exception to RuggedException
            DummyLocalizable message = new DummyLocalizable(e.getLocalizedMessage());
            final RuggedException re = new RuggedException(message, this.geoidFilePath, e.getLocalizedMessage());
            re.initCause(e);
            throw re;
        }
    }

    /**
     * Read the Geoid dataset (initialize elevationsWholeEarth) and initialize Geoid sub-tiles (according to the Geoid data read)
     * @throws S2GeolibException
     */
    private void initGeoidSubTiles() throws RuggedException {

        // Open the Geoid dataset
        // ======================
        Dataset geoidDataset = openGeoid();

        // Analyze and read the raster: initialize this.elevationsWholeEarth
        // ===========================
        double[] geoTransformInfo = readRaster(geoidDataset);

        // Define the 8 Geoid sub-tiles
        // ============================
        defineGeoidSubTiles(this.rasterLonSize, this.rasterLatSize, geoTransformInfo);

    }

    /**
     * Analyse and read the Geoid raster. Initialize elevationsWholeEarth (fill in with GDAL convention : [0][0] = North West)
     * @param geoidDataset the geoid Dataset (opened by GDAL)
     * @return the geoTransformInfo (read from the geoidDataset)
     * @throws S2GeolibException
     */
    private double[] readRaster(Dataset geoidDataset) throws RuggedException {

        // Analyse the Geoid dataset
        // =========================
        this.rasterLonSize = geoidDataset.GetRasterXSize();
        this.rasterLatSize = geoidDataset.GetRasterYSize();

        // Get the Geo transform info: lon origin (deg), lon step (deg), 0, lat origin (deg), 0, lat step (deg)
        double[] geoTransformInfo = geoidDataset.GetGeoTransform();
        if (geoTransformInfo.length < 6) { // should not occur; if a problem occurs : the default transform value is set to (0,1,0,0,0,1)
            String str = "GDALGeoTransform does not contains 6 elements";
            DummyLocalizable message = new DummyLocalizable(str);
            throw new RuggedException(message);
        }

        // Read the raster
        // ===============
        // To store the raster data
        double[] data = new double[rasterLatSize * rasterLonSize];
        // Get the raster
        Band band = geoidDataset.GetRasterBand(1);
        // We read all the raster once for all
        band.ReadRaster(0, 0, rasterLonSize, rasterLatSize, data);
        // Get the no data value from the raster
        Double[] noDataValue = new Double[1];
        band.GetNoDataValue(noDataValue);
        // test if the no data value exists
        Boolean noDataValueExist = false;
        if (noDataValue[0] != null ) noDataValueExist = true;

        // Read the full raster all in once ...
        this.elevationsWholeEarth = new double[rasterLatSize][rasterLonSize];

        // from bottom left to upper right corner (GDAL convention)
        for (int iLat = rasterLatSize - 1; iLat >= 0; iLat--) {
            for (int jLon = 0; jLon < rasterLonSize; jLon++) {
                // if jlon = 0 and ilat = rasterLatSize (721) : crash
                // if iLat = 720 and jlon = rasterLonSize (1440) : crash

                // Read the data with GDAL convention : iLat = 0 at North; iLat = rasterLatSize - 1 at South; jLon = 0 at West; jLon = rasterLonSize - 1 at East
                double elevationOverEllipsoid = data[iLat * rasterLonSize + jLon];
                if (noDataValueExist && (elevationOverEllipsoid == noDataValue[0])){
                    elevationOverEllipsoid = 0.0;
                }
                // Set elevation for current latitude and longitude
                // IMPORTANT : keep the GDAL convention [0][0] = North West
                this.elevationsWholeEarth[iLat][jLon] = elevationOverEllipsoid;

            } // end loop jLon
        } // end loop iLat


        // Close Dataset and Band once the raster is read ...
        band.delete();
        band = null;
        geoidDataset.delete();
        geoidDataset = null;

        // Return the Geo Transform info
        return geoTransformInfo;
    }

    /**
     * Split of the Earth in 8 tiles : 2 bands in latitude and 4 bands in longitude.
     * Each tile must have 1 pixel overlapping with the neighbor.
     * @param rasterLonSize size of the raster in longitude
     * @param rasterLatSize size of the raster in latitude
     * @param geoTransformInfo the geoTransformInfo read from the geoid Dataset: lon origin (deg), lon step (deg), 0, lat origin (deg), 0, lat step (deg)
     * @throws S2GeolibException
     */
    private void defineGeoidSubTiles(int rasterLonSize, int rasterLatSize, double[] geoTransformInfo) throws RuggedException {

        GeoTransformTools gtTools = new GeoTransformTools(geoTransformInfo, rasterLonSize, rasterLatSize);
        this.lonStepInDeg = FastMath.abs(geoTransformInfo[1]);
        this.latStepInDeg = FastMath.abs(geoTransformInfo[5]);

        //          Size is 1440, 721
        //          Origin = (-180.125000000000000,90.125000000000000)
        //                Pixel Size = (0.250000000000000,-0.250000000000000)
        //                Corner Coordinates:
        //                Upper Left  (-180.1250000,  90.1250000) (180d 7'30.00"W, 90d 7'30.00"N)
        //                Lower Left  (-180.1250000, -90.1250000) (180d 7'30.00"W, 90d 7'30.00"S)
        //                Upper Right ( 179.8750000,  90.1250000) (179d52'30.00"E, 90d 7'30.00"N)
        //                Lower Right ( 179.8750000, -90.1250000) (179d52'30.00"E, 90d 7'30.00"S)
        //                Center      (  -0.1250000,   0.0000000) (  0d 7'30.00"W,  0d 0' 0.01"N)

        // For Latitude : step is negative ...
        // ============
        // Latitude band North: from North pole to equator
        // Latitude band South: from equator to South pole
        // ----------------------------------------------------------------------------------------------------
        // Raster size from 90.125 (North) to -90.125 (South) = 1 pixel of overlapping for latitude in the data
        // ----------------------------------------------------------------------------------------------------

        // Compute band latitudes from line numbers (latitude corresponds to the upperleft corner of the DEM cell)
        int nbRowsLatBand = rasterLatSize/2; // 2 bands in latitude; Must be < rasterLatSize !!

        // North band :
        int lineLatMaxNorth = 0;
        int lineLatMinNorth = nbRowsLatBand + 1; // 1 line of overlapping with South band
        // Not use in fact : double latMaxNorth = gtTools.getYFromPixelLine(0, lineLatMaxNorth);
        double latMinNorth = gtTools.getYFromPixelLine(0, lineLatMinNorth); // 1 line of overlapping with South band
        int latSizeNorth = 1 + (lineLatMinNorth - lineLatMaxNorth);

        // South band:
        int lineLatMaxSouth = nbRowsLatBand; // Line of overlapping taking into account in lineLatMinNorth
        int lineLatMinSouth = rasterLatSize - 1; // When reading the data : cannot go up to rasterLatSize
        double latMaxSouth = gtTools.getYFromPixelLine(0, lineLatMaxSouth); // Line of overlapping taking into account in latMinNorth
        double latMinSouth = gtTools.getYFromPixelLine(0, lineLatMinSouth);
        int latSizeSouth= 1 + (lineLatMinSouth - lineLatMaxSouth);

        // For longitude : step is positive
        // =============
        // Longitude band FarWest : from West anti-meridian (-180) to -90 1439-1080
        // Longitude band West : from -90 to Greenwich meridian (0)
        // Longitude band East : from Greenwich meridian (0) to +90
        // Longitude band FarEast: from +90 to East anti-meridian (+180)
        // ----------------------------------------------------------------------------------------------------
        // Raster size from -180.125 (West) to 179.875 (East) = No pixel of overlapping for longitude in the data
        // ----------------------------------------------------------------------------------------------------

        // Compute band longitudes from line numbers (longitude corresponds to the upperleft corner of the DEM cell)
        int nbColsLonBand = rasterLonSize / 4; // 4 bands in longitude; Must be < rasterLonSize !!

        // FarWest band :
        int pxLonMinFarWest = 0;
        int pxLonMaxFarWest = nbColsLonBand + 1; // 1 column of overlapping with West band
        // The Far West min longitude = -180.125
        double lonMinFarWest = gtTools.getXFromPixelLine(pxLonMinFarWest, 0);
        double lonMaxFarWest = gtTools.getXFromPixelLine(pxLonMaxFarWest, 0);
        int lonSizeFarWest = 1 + (pxLonMaxFarWest - pxLonMinFarWest);

        // West band :
        int pxLonMinWest = nbColsLonBand;        // Line of overlapping taken into account in pxLonMaxFarWest
        int pxLonMaxWest = 2*nbColsLonBand + 1;  // 1 column of overlapping with East band
        double lonMinWest = gtTools.getXFromPixelLine(pxLonMinWest, 0);
        double lonMaxWest = gtTools.getXFromPixelLine(pxLonMaxWest, 0);
        int lonSizeWest = 1 + (pxLonMaxWest - pxLonMinWest);

        // East band :
        int pxLonMinEast = 2*nbColsLonBand;     // Line of overlapping taken into account in pxLonMaxWest
        int pxLonMaxEast = 3*nbColsLonBand + 1; // 1 column of overlapping with FarEast band
        double lonMinEast = gtTools.getXFromPixelLine(pxLonMinEast, 0);
        double lonMaxEast = gtTools.getXFromPixelLine(pxLonMaxEast, 0);
        int lonSizeEast = 1 + (pxLonMaxEast - pxLonMinEast);

        // FarEast band:
        int pxLonMinFarEast = 3*nbColsLonBand;    // Line of overlapping taken into account in pxLonMaxEast
        int pxLonMaxFarEast = rasterLonSize - 1;  // When reading the data : cannot go up to rasterLonSize
        double lonMinFarEast = gtTools.getXFromPixelLine(pxLonMinFarEast, 0);
        // The Far East max longitude = 179.875 ; no overlapping with the Far West band that starts at -180.125
        // When loading the elevations: this band will need the first column of the FarWest band
        // Not use in fact : double lonMaxFarEast = gtTools.getXFromPixelLine(pxLonMaxFarEast, 0);
        // We need to overlap the tile at anti-meridian (between +180 deg and -180 deg).
        // We add 1 pixel to define the correct tile size ...
        // The pxLonMaxFarEast doesn't reach the extra pixel (that will be filled in by the FarWest corresponding tile)
        int lonSizeFarEast = 1 + ( 1 + (pxLonMaxFarEast - pxLonMinFarEast));

        // Initialize the tiles definition
        // ===============================
        this.tile_NORTH_FARWEST = new GeoidTileDefinition(EarthQuadrant.NORTH_FARWEST,
                                                          latMinNorth, lonMinFarWest, latSizeNorth, lonSizeFarWest,
                                                          lineLatMinNorth, lineLatMaxNorth, pxLonMinFarWest, pxLonMaxFarWest);
        this.tile_NORTH_WEST    = new GeoidTileDefinition(EarthQuadrant.NORTH_WEST,
                                                          latMinNorth, lonMinWest,    latSizeNorth, lonSizeWest,
                                                          lineLatMinNorth, lineLatMaxNorth, pxLonMinWest, pxLonMaxWest);
        this.tile_NORTH_EAST    = new GeoidTileDefinition(EarthQuadrant.NORTH_EAST,
                                                          latMinNorth, lonMinEast,    latSizeNorth, lonSizeEast,
                                                          lineLatMinNorth, lineLatMaxNorth, pxLonMinEast, pxLonMaxEast);
        this.tile_NORTH_FAREAST = new GeoidTileDefinition(EarthQuadrant.NORTH_FAREAST,
                                                          latMinNorth, lonMinFarEast, latSizeNorth, lonSizeFarEast,
                                                          lineLatMinNorth, lineLatMaxNorth, pxLonMinFarEast, pxLonMaxFarEast);

        this.tile_SOUTH_FARWEST = new GeoidTileDefinition(EarthQuadrant.SOUTH_FARWEST,
                                                          latMinSouth, lonMinFarWest, latSizeSouth, lonSizeFarWest,
                                                          lineLatMinSouth, lineLatMaxSouth, pxLonMinFarWest, pxLonMaxFarWest);
        this.tile_SOUTH_WEST    = new GeoidTileDefinition(EarthQuadrant.SOUTH_WEST,
                                                          latMinSouth, lonMinWest,    latSizeSouth, lonSizeWest,
                                                          lineLatMinSouth, lineLatMaxSouth, pxLonMinWest, pxLonMaxWest);
        this.tile_SOUTH_EAST    = new GeoidTileDefinition(EarthQuadrant.SOUTH_EAST,
                                                          latMinSouth, lonMinEast,    latSizeSouth, lonSizeEast,
                                                          lineLatMinSouth, lineLatMaxSouth, pxLonMinEast, pxLonMaxEast);
        this.tile_SOUTH_FAREAST = new GeoidTileDefinition(EarthQuadrant.SOUTH_FAREAST,
                                                          latMinSouth, lonMinFarEast, latSizeSouth, lonSizeFarEast,
                                                          lineLatMinSouth, lineLatMaxSouth, pxLonMinFarEast, pxLonMaxFarEast);

        // Initialize the Earth quadrant limits
        // ====================================
        this.earthQuadrantsLimits = new EarthQuadrantsLimits(latMinNorth, latMaxSouth, lonMaxFarWest, lonMaxWest, lonMaxEast);
    }

    /**
     * Identify the tile according to the Earth quadrant
     * @param currentQuadrant Earth quadrant where the (latitude, longitude) belongs to
     * @return the Geoid tile definition according to the quadrant
     */
    private GeoidTileDefinition identifyGeoidTile(EarthQuadrant currentQuadrant) {

        // Search which tile the (latitude, longitude) belongs to
        // ======================================================
        GeoidTileDefinition foundGeoidTile = null;

        switch (currentQuadrant) {
            case NORTH_FARWEST:
                foundGeoidTile = this.tile_NORTH_FARWEST;
                break;

            case NORTH_WEST:
                foundGeoidTile = this.tile_NORTH_WEST;
                break;

            case NORTH_EAST:
                foundGeoidTile = this.tile_NORTH_EAST;
                break;

            case NORTH_FAREAST:
                foundGeoidTile = this.tile_NORTH_FAREAST;
                break;

            case SOUTH_FARWEST:
                foundGeoidTile = this.tile_SOUTH_FARWEST;
                break;

            case SOUTH_WEST:
                foundGeoidTile = this.tile_SOUTH_WEST;
                break;

            case SOUTH_EAST:
                foundGeoidTile = this.tile_SOUTH_EAST;
                break;

            case SOUTH_FAREAST:
                foundGeoidTile = this.tile_SOUTH_FAREAST;
                break;

            default:
                break;
        }

        return foundGeoidTile;
    }

    /**
     * Fill in the tile
     * @param foundGeoidTile Geoid sub-tile found where the (latitude, longitude) belongs to
     * @param tile tile to update
     * @return the last longitude index filled in for the tile (needed for FarEast tile to complete the last column)
     * @throws RuggedException
     */
    private int fillTheTile(GeoidTileDefinition foundGeoidTile, UpdatableTile tile) throws RuggedException {

        // Init Tile Geometry for Rugged
        // -----------------------------
        // For Rugged, minimum latitude and longitude correspond to the CENTRE of the farthest SouthWest cell of the DEM.
        // For GDAL, the latitude and longitude correspond to the UPPER LEFT coordinate of the cell of the DEM.
        double latMinInRad = FastMath.toRadians(-0.5*this.latStepInDeg + foundGeoidTile.getLatMinDeg());
        double lonMinInRad = FastMath.toRadians(0.5*this.lonStepInDeg + foundGeoidTile.getLonMinDeg());

        int latRows = foundGeoidTile.getLatRows();
        int lonColumns = foundGeoidTile.getLonColumns();
        double latStepInRad = FastMath.toRadians(this.latStepInDeg);
        double lonStepInRad = FastMath.toRadians(this.lonStepInDeg);

        tile.setGeometry(latMinInRad, lonMinInRad, latStepInRad, lonStepInRad, latRows, lonColumns);

        // Fill in the tile
        // ----------------
        // For raster : latitude line = 0 at North and 720 at South ;
        //       Here : lineLatMax = North line number and lineLatMin = South line number ==> inversion of loop !
        int iSouthTile = foundGeoidTile.getLineLatMin();
        int iNorthTile = foundGeoidTile.getLineLatMax();

        int jLonTile = Integer.MIN_VALUE;
        for (int iLat = iNorthTile; iLat <= iSouthTile; iLat++) {

            // For Rugged : line = 0 at South (reverse from the raster numbering)
            int iLatTile = iSouthTile - iLat; // Index for latitude to fill in the tile

            // Reset the longitude index ...
            jLonTile = -1;

            // For FarEast tiles : the tile size has 1 extra column (lonColumns value) than the range PxLonMin to PxLonMax ...
            //                     The FarEast tiles will be completed in fillTheLastColumnTile
            for (int jLon = foundGeoidTile.getPxLonMin(); jLon <= foundGeoidTile.getPxLonMax(); jLon++) {
                jLonTile++; // Index for longitude to fill in the tile
                double elevationOverEllipsoid = this.elevationsWholeEarth[iLat][jLon];
                // Fill the tile at the proper index ...
                tile.setElevation(iLatTile, jLonTile, elevationOverEllipsoid);
            } // end loop jLon
        } // end loop iLat

        // Return the last longitude index of the tile
        return jLonTile;
    }

    /**
     * Fill the last column of the FarEast tiles with the first column of the corresponding FarWest tile
     * @param farWestGeoidTile FarWest Geoid sub-tile where to read the first column
     * @param tile FarEast tile with the last column to update (because it is a FarEast tile)
     * @param indexLonMax index max in longitude found for the FarEast tile
     * @throws RuggedException
     */
    private void fillTheLastColumnTile(GeoidTileDefinition farWestGeoidTile, UpdatableTile tile, int indexLonMax) throws RuggedException {

        // For raster : latitude line = 0 at North and 720 at South ;
        // Here : lineLatMax = North line number and lineLatMin = South line number ==> inversion of loop !
        int iSouthTile = farWestGeoidTile.getLineLatMin();
        int iNorthTile = farWestGeoidTile.getLineLatMax();

        // compute the last column from the index max in longitude found for the FarEast tile
        int jLonLast = indexLonMax + 1;

        for (int iLat = iNorthTile; iLat <= iSouthTile; iLat++) {
            // Get the first column elevation from the FarWestGeoidTile
            double elevationOverEllipsoid = this.elevationsWholeEarth[iLat][0];

            // For Rugged : line = 0 at South (reverse from the raster numbering)
            int iLatTile = iSouthTile - iLat; // Index for latitude to fill in the tile
            // Fill the tile at the proper index ...
            tile.setElevation(iLatTile, jLonLast, elevationOverEllipsoid);
        } // end loop iLat
    }

    /**
     * Search the Earth quadrant where the (latitude, longitude) belongs to
     * @param latitude normalized latitude between +/- 90 deg (deg)
     * @param longitude normalized longitude between +/- 180 deg (deg)
     * @return the EarthQuadrant defines by the latitude and the longitude
     */
    private EarthQuadrant getQuadrant(double latitude, double longitude) {

        EarthHemisphere zoneLat= null;
        if (latitude < this.earthQuadrantsLimits.getLimitSouthNorth()) { // South zone
            //         NOT compatible with   Rugged SimpleTile.getLocation search for HAS_INTERPOLATION_NEIGHBORS:
            //         if (latitude < this.earthQuadrantsLimits.getLatMaxSouthDeg()) { // South zone

            zoneLat = EarthHemisphere.SOUTH;
        } else { // North zone
            zoneLat = EarthHemisphere.NORTH;
        }

        EarthHemisphere zoneLon = null;
        if (longitude < this.earthQuadrantsLimits.getLonMaxWestDeg()){ // West from Greenwich
            if (longitude < this.earthQuadrantsLimits.getLonMaxFarWestDeg()){ // FarWest zone
                zoneLon = EarthHemisphere.FARWEST;
            } else { // West zone
                zoneLon = EarthHemisphere.WEST;
            }
        } else { // East from Greenwich
            if (longitude < this.earthQuadrantsLimits.getLonMaxEastDeg()) { // East zone
                zoneLon = EarthHemisphere.EAST;
            } else { // FarEast zone
                zoneLon = EarthHemisphere.FAREAST;
            }
        }
        return getQuadrant(zoneLat, zoneLon);
    }

    /**
     * Search the Earth quadrant where the latitude EarthHemisphere and longitude EarthHemisphere belongs to
     * @param zoneLat EarthHemisphere for latitude (North or South)
     * @param zoneLon EarthHemisphere for longitude (FarWest, West, East or FarEast)
     * @return the EarthQuadrant defines by the latitude EarthHemisphere and longitude EarthHemisphere
     */
    private EarthQuadrant getQuadrant(EarthHemisphere zoneLat, EarthHemisphere zoneLon){
        EarthQuadrant quadrant = null;
        if (zoneLat == EarthHemisphere.SOUTH){
            if (zoneLon == EarthHemisphere.FARWEST) quadrant = EarthQuadrant.SOUTH_FARWEST;
            if (zoneLon == EarthHemisphere.WEST) quadrant = EarthQuadrant.SOUTH_WEST;
            if (zoneLon == EarthHemisphere.EAST) quadrant = EarthQuadrant.SOUTH_EAST;
            if (zoneLon == EarthHemisphere.FAREAST) quadrant = EarthQuadrant.SOUTH_FAREAST;
        } else { // zoneLat = Zone.NORTH
            if (zoneLon == EarthHemisphere.FARWEST) quadrant = EarthQuadrant.NORTH_FARWEST;
            if (zoneLon == EarthHemisphere.WEST) quadrant = EarthQuadrant.NORTH_WEST;
            if (zoneLon == EarthHemisphere.EAST) quadrant = EarthQuadrant.NORTH_EAST;
            if (zoneLon == EarthHemisphere.FAREAST) quadrant = EarthQuadrant.NORTH_FAREAST;
        }
        return quadrant;
    }

    /**
     * Open the Geoid dataset
     * @throws S2GeolibException
     */
    public Dataset openGeoid() throws RuggedException {
        Dataset geoidDataset = gdal.Open(geoidFilePath, gdalconst.GA_ReadOnly);

        if (geoidDataset == null) {
            String str = geoidFilePath + "cannot be read";
            DummyLocalizable message = new DummyLocalizable(str);
            throw new RuggedException(message);
        }
        return geoidDataset;
    }

    //   protected double getCoordInsideFootprint(double coord, double minCoord, double maxCoord, double maxBorderValue) {
    //      double returned = coord;
    //      if (coord > maxCoord) {
    //         // Take distance from border
    //         double distanceFromBorder = maxBorderValue - coord;
    //         // Get new coord from other border
    //         returned = -maxBorderValue - distanceFromBorder;
    //      } else if (coord < minCoord) {
    //         // coord < minCoord
    //         double distanceFromBorder = -maxBorderValue - coord;
    //         returned = maxBorderValue + distanceFromBorder;
    //      }
    //      return returned;
    //   }

    //    /**
    //     * Return Geoid longitude step ijLatn degrees
    //     * @param latitude (rad)
    //     * @param longitude (rad)
    //     * @return
    //     * @throws RuggedException
    //     */
    //    public double getGeoidStep(double latitude, double longitude) throws S2GeolibException {
    //
    //        double lonStep = Double.NaN;
    //
    //        if (geoidDataset == null) {
    //            openGeoid();
    //        }
    //
    //        double[] geoTransformInfo = geoidDataset.GetGeoTransform();
    //
    //        if (geoTransformInfo.length < 6) {
    //            String[] strParams = new String[]{geoidFilePath, "GDAL geoTransformInfo < 6"};
    //            String message = S2GeolibResourceBundle.getString(S2GeolibResourceBundle.ERROR_PARSING_FILE, strParams);
    //            throw new S2GeolibException(message);
    //        } else {
    //            lonStep = FastMath.abs(geoTransformInfo[1]);
    //        }
    //        return lonStep;
    //    }

    //    /**
    //     * Return the full Geoid Dataset
    //     * @return the geoid dataset
    //     * @throws RuggedException
    //     */
    //    protected Dataset getGeoidDataset() throws RuggedException {
    //
    //        if (this.geoidDataset == null) {
    //            try {
    //                openGeoid();
    //            } catch (S2GeolibException e) {
    //                DummyLocalizable message = new DummyLocalizable(e.getMessage());
    //                RuggedException re = new RuggedException(message, e);
    //                throw re;
    //            }
    //        }
    //        return this.geoidDataset;
    //    }

    //    /**
    //     * @return the geoidFilePath
    //     */
    //    public String getGeoidFilePath() {
    //        return geoidFilePath;
    //    }
}
