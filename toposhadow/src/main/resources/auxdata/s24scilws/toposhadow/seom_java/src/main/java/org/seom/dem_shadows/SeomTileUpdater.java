package org.seom.dem_shadows;

import org.orekit.rugged.raster.SimpleTile;
import org.orekit.rugged.raster.SimpleTileFactory;
import org.orekit.rugged.raster.TileUpdater;
import org.orekit.rugged.raster.TilesCache;

public abstract class SeomTileUpdater implements TileUpdater {

    /**
     * Tiles cache used to improve perf
     */
    protected TilesCache<SimpleTile> tilesCache = null;

    /**
     * Constructor
     */
    public SeomTileUpdater() {
        SimpleTileFactory tileFactory = new SimpleTileFactory();
//        tilesCache = new TilesCache<>(tileFactory, this, S2GeoConstants.NB_CACHE_TILE);
        tilesCache = new TilesCache<>(tileFactory, this, 8); // TODO
    }

//    /**
//     * Get elevation at given lat lon (lat/lon is expected to be in same srs than DEM)
//     * @param latitude (degree)
//     * @param longitude (degree)
//     * @return
//     * @throws S2GeoException
//     */
//    public double getElevation(double latitude, double longitude) {
//        
//    	double elevation = Double.NaN;
//        elevation = getElevationRadian(latitude, longitude);
//
//        return elevation;
//    }
//
//    /**
//     * Get elevation at given lat lon (lat/lon is expected to be in same srs than DEM)
//     * @param latitude (radian)
//     * @param longitude (radian)
//     * @return
//     * @throws S2GeoException
//     */
//    public double getElevationRadian(double latitude, double longitude) {
//        double elevation = Double.NaN;
//        try {
//            SimpleTile tile = tilesCache.getTile(latitude, longitude);
//            elevation = tile.interpolateElevation(latitude, longitude);
//        } catch (Exception e) {
//        }
//        return elevation;
//    }
}
