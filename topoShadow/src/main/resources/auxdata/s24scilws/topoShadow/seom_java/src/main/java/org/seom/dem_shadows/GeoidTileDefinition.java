package org.seom.dem_shadows;



/**
 * Definition of one tile of Geoid
 * @author Guylaine Prat
 *
 */
public class GeoidTileDefinition {

    /**
     * Earth quadrant where the Geoid tile belongs to
     */
    private EarthQuadrant earthQuadrant = null;

    /**
     * Latitude min of the tile (deg)
     * TBN: this latitude corresponds to the upper left corner of the DEM cell (defined by lineLatMin and pxLonMin)
     */
    private double latMinDeg = Double.NaN;

    /**
     * Longitude min of the tile (deg)
     * TBN: this longitude cooresponds to the upper left corner of the DEM cell (defined by pxLonMin and lineLatMin)
     */
    private double lonMinDeg = Double.NaN;

    /**
     * Number of rows in latitude of the tile
     */
    private int latRows = Integer.MAX_VALUE;

    /**
     * Number of columns in longitude of the tile
     */
    private int lonColumns = Integer.MAX_VALUE;

    /**
     * Associated line number to the latMinDeg
     */
    private int lineLatMin = Integer.MAX_VALUE;

    /**
     * Associated line number to the maximum latitude of the tile
     */
    private int lineLatMax = Integer.MIN_VALUE;

    /**
     * Associated column number to the lonMinDeg
     */
    private int pxLonMin = Integer.MAX_VALUE;

    /**
     * Associated column number to the maximum longitude of the tile
     */
    private int pxLonMax =  Integer.MIN_VALUE;


    /**
     * @param earthQuadrant
     * @param latMinDeg
     * @param lonMinDeg
     * @param latRows
     * @param lonColumns
     * @param lineLatMin
     * @param lineLatMax
     * @param pxLonMin
     * @param pxLonMax
     */
    public GeoidTileDefinition(EarthQuadrant earthQuadrant, double latMinDeg, double lonMinDeg,
                               int latRows, int lonColumns,
                               int lineLatMin, int lineLatMax,
                               int pxLonMin, int pxLonMax) {
        this.earthQuadrant = earthQuadrant;
        this.latMinDeg = latMinDeg;
        this.lonMinDeg = lonMinDeg;
        this.latRows = latRows;
        this.lonColumns = lonColumns;
        this.lineLatMin = lineLatMin;
        this.lineLatMax = lineLatMax;
        this.pxLonMin = pxLonMin;
        this.pxLonMax = pxLonMax;
    }

    /**
     * Get the Earth quadrant the tile belongs to
     * @return Earth quadrant the tile belongs to
     */
    public EarthQuadrant getEarthQuadrant() {
        return earthQuadrant;
    }

    /**
     * Get the latitude min of the tile (deg)
     * TBN: this latitude corresponds to the upper left corner of the DEM cell (defined by lineLatMin and pxLonMin)
     * @return latitude min of the tile (deg)
     */
    public double getLatMinDeg() {
        return latMinDeg;
    }

    /**
     * Get the longitude min of the tile (deg)
     * TBN: this longitude corresponds to the upper left corner of the DEM cell (defined by pxLonMin and lineLatMin)
     * @return longitude min of the tile (deg)
     */
    public double getLonMinDeg() {
        return lonMinDeg;
    }

    /**
     * Get the number of rows in latitude of the tile
     * @return  number of rows in latitude of the tile
     */
    public int getLatRows() {
        return latRows;
    }

    /**
     * Get the number of columns in longitude of the tile
     * @return number of columns in longitude of the tile
     */
    public int getLonColumns() {
        return lonColumns;
    }

    /**
     * Get the associated line number to the latMinDeg
     * @return associated line number to the latMinDeg
     */
    public int getLineLatMin() {
        return lineLatMin;
    }

    /**
     * Get the associated line number to the maximum latitude of the tile
     * @return associated line number to the maximum latitude of the tile
     */
    public int getLineLatMax() {
        return lineLatMax;
    }

    /**
     * Get the associated column number to the lonMinDeg
     * @return associated column number to the lonMinDeg
     */
    public int getPxLonMin() {
        return pxLonMin;
    }

    /**
     * Get the associated column number to the maximum longitude of the tile
     * @return associated column number to the maximum longitude of the tile
     */
    public int getPxLonMax() {
        return pxLonMax;
    }
}
