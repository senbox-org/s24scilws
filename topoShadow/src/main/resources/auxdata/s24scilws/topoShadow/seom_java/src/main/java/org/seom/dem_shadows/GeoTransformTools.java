package org.seom.dem_shadows;



/**
 * Tool using GeoTransform info returned from Gdal
 * @author Guylaine Prat
 */
public class GeoTransformTools {

    /**
     * GeoTransform info from Gdal
     */
    private double[] gtInfo;
    /**
     * x size (longitude in deg)
     */
    private double xSize;

    /**
     * y size (latitude in deg)
     */
    private double ySize;

    /**
     * @param gtInfo geo transform info from gdal
     * @param xSize x size (longitude in deg)
     * @param ySize y size (latitude in deg)
     */
    public GeoTransformTools(double[] gtInfo, double xSize, double ySize) {
        // We suppose gtInfo[5] and gtInfo[1] different from 0
        // no tests are performed about gtInfo values ... as org.gdal.gdal.Dataset.GetGeoTransform will give a default transform value set to (0,1,0,0,0,1) if problems
        this.gtInfo = gtInfo;
        this.xSize = xSize;
        this.ySize = ySize;
    }

    /**
     * Get X (longitude) coord from pixel and line
     * @param pixelNum pixel number
     * @param lineNum line number
     * @return X coord (longitude in degrees)  from pixel and line. For Gdal = Upper left corner of the cell.
     */
    public double getXFromPixelLine(double pixelNum, double lineNum) {
        return gtInfo[0] + gtInfo[1] * pixelNum + gtInfo[2] * lineNum;
    }

    /**
     * Get Y (latitude) coord from pixel and line
     * @param pixelNum pixel number
     * @param lineNum line number
     * @return Y coord (latitude in degrees) from pixel and line. For Gdal = Upper left corner of the cell.
     */
    public double getYFromPixelLine(double pixelNum, double lineNum) {
        return gtInfo[3] + gtInfo[4] * pixelNum + gtInfo[5] * lineNum;
    }

    /**
     * Get pixel number from latitude and longitude
     * @param lat latitude (deg)
     * @param lon longitude (deg)
     * @return pixel number. . For Gdal = Upper left corner of the cell.
     */
    public double getPixelNumFromLatLon(double lat, double lon) {
        // We suppose gtInfo[5] and gtInfo[1] different from 0
        double tmp = gtInfo[2] / gtInfo[5];
        double value = (lon - tmp * lat - gtInfo[0] + tmp * gtInfo[3]) / (gtInfo[1] - tmp * gtInfo[4]);
        if (value >= xSize) { // if value = xSize, don't convert value to get the inverse value from getLatMax and getLonMax
            value = value % xSize;
        }
        return value;

        //return (lon - gtInfo[0]) / gtInfo[1];
    }

    /**
     * Get line number from latitude and longitude
     * @param lat latitude (deg)
     * @param lon longitude (deg)
     * @return line number. For Gdal = Upper left corner of the cell.
     */
    public double getLineNumFromLatLon(double lat, double lon) {
        // We suppose gtInfo[1] and gtInfo[5] different from 0
        double tmp = gtInfo[4] / gtInfo[1];
        double value = (lat - tmp * lon - gtInfo[3] + tmp * gtInfo[0]) / (gtInfo[5] - tmp * gtInfo[2]);
        if (value >= ySize) { // if value = xSize, don't convert value to get the inverse value from getLatMax and getLonMax
            value = value % ySize;
        }
        return value;
        //return (lat - gtInfo[3]) / gtInfo[5];
    }

    //    /**
    //     * Convert the given degrees value in meters
    //     * @param distanceDeg
    //     * @return the given degrees value in meters
    //     */
    //    public static double computeDistanceDeg(double lonDeg, double latDeg, double distanceDeg) {
    //        double lonDeg1 = lonDeg;
    //        double latDeg1 = latDeg;
    //        double lonDeg2 = lonDeg + distanceDeg;
    //        double latDeg2 = latDeg;
    //        return computeDistance(lonDeg1, latDeg1, lonDeg2, latDeg2);
    //    }
    //
    //    /**
    //     * Compute distance between point (lonDeg1, latDeg1) and point (lonDeg2, latDeg2) in meters
    //     * @param lonDeg1
    //     * @param latDeg1
    //     * @param lonDeg2
    //     * @param latDeg2
    //     * @return distance between point (lonDeg1, latDeg1) and point (lonDeg2, latDeg2)
    //     */
    //    public static double computeDistance(double lonDeg1, double latDeg1, double lonDeg2, double latDeg2) {
    //        double returned;
    //        double xRad1 = FastMath.toRadians(lonDeg1);
    //        double xRad2 = FastMath.toRadians(lonDeg2);
    //        double yRad1 = FastMath.toRadians(latDeg1);
    //        double yRad2 = FastMath.toRadians(latDeg2);
    //        returned = computeDistanceRad(xRad1, yRad1, xRad2, yRad2);
    //        return returned;
    //    }
    //
    //    /**
    //     * distance in meters between point (xRad1, yRad1) and point (xRad2, yRad2)
    //     * @param xRad1
    //     * @param xRad2
    //     * @param yRad1
    //     * @param yRad2
    //     * @return
    //     */
    //    private static double computeDistanceRad(double xRad1, double yRad1, double xRad2, double yRad2) {
    //        double returned;
    //        double sinValue = FastMath.sin(xRad1) * FastMath.sin(xRad2);
    //        double deltaLambda = FastMath.abs(yRad1 - yRad2);
    //        double cosValue = FastMath.cos(xRad1) * FastMath.cos(xRad2) * FastMath.cos(deltaLambda);
    //        double deltaPhy = FastMath.acos(sinValue + cosValue);
    //        if (Double.isNaN(deltaPhy)) {
    //            deltaPhy = 0d;
    //        }
    //        returned = S2GeolibConstants.EARTH_RADIUS / 100 * deltaPhy;
    //        return returned;
    //    }

    public double getLonMin() {
        return gtInfo[0];
    }

    public double getLatMin() {
        return gtInfo[3] + xSize * gtInfo[4] + gtInfo[5] * ySize;
    }

    public double getLonMax() {
        return gtInfo[0] + xSize * gtInfo[1] + gtInfo[2] * ySize;
    }

    public double getLatMax() {
        return gtInfo[3];
    }

    //    /**
    //     * Create a Geometry object with a margin. This geometry will be used to filter detector mask
    //     * @param maxY
    //     * @param minX
    //     * @param maxX
    //     * @param minY
    //     * @param yStep
    //     * @param xStep
    //     * @return
    //     */
    //    public static Geometry createTileGeometry(double maxY, double minX, double maxX, double minY, double yStep, double xStep,
    //           CoordinateTransformation transformer) {
    //        double xMin = minX - 15 * xStep;
    //        double xMax = maxX + 15 * xStep;
    //        double yMin = minY - 15 * yStep;
    //        double yMax = maxY + 15 * yStep;
    //        // Load all detector geometry to reduce inverse location time
    //
    //        // Create filter with tile geometry
    //        Geometry filterPoPolygon = new Geometry(ogr.wkbPolygon);
    //
    //        Geometry filterRingGeom = new Geometry(ogrConstants.wkbLinearRing);
    //
    //        // Add coords point by point.
    //        // We must call TransformPoint on every corner because of rotation between SRS
    //        double[] transformed = { xMin, yMin };
    //        if (transformer != null) {
    //            transformed = transformer.TransformPoint(transformed[0], transformed[1]);
    //        }
    //        // Save first point coords to close the ring
    //        double[] firstPoint = transformed;
    //        filterRingGeom.AddPoint(firstPoint[0], firstPoint[1]);
    //
    //        transformed = new double[] { xMin, yMax };
    //        if (transformer != null) {
    //            transformed = transformer.TransformPoint(transformed[0], transformed[1]);
    //        }
    //        filterRingGeom.AddPoint(transformed[0], transformed[1]);
    //
    //        transformed = new double[] { xMax, yMax };
    //        if (transformer != null) {
    //            transformed = transformer.TransformPoint(transformed[0], transformed[1]);
    //        }
    //        filterRingGeom.AddPoint(transformed[0], transformed[1]);
    //
    //        transformed = new double[] { xMax, yMin };
    //        if (transformer != null) {
    //            transformed = transformer.TransformPoint(transformed[0], transformed[1]);
    //        }
    //        filterRingGeom.AddPoint(transformed[0], transformed[1]);
    //        // Close the ring with first point coords
    //        filterRingGeom.AddPoint(firstPoint[0], firstPoint[1]);
    //        filterPoPolygon.AddGeometry(filterRingGeom);
    //        return filterPoPolygon;
    //    }

}

