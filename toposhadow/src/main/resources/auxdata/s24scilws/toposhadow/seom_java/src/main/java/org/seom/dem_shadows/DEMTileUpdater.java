package org.seom.dem_shadows;


import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.rugged.errors.RuggedException;
import org.orekit.rugged.raster.SimpleTile;
import org.orekit.rugged.raster.SimpleTileFactory;
import org.orekit.rugged.raster.TileUpdater;
import org.orekit.rugged.raster.TilesCache;


public abstract class DEMTileUpdater implements TileUpdater {

    /**
     * Tiles cache used to improve perf
     */
    private TilesCache<SimpleTile> tilesCache = null;

    /**
     * Constructor
     */
    public DEMTileUpdater() {
        SimpleTileFactory tileFactory = new SimpleTileFactory();
        tilesCache = new TilesCache<>(tileFactory, this, Constants.NB_CACHE_TILE);
    }

    /**
     * Get elevation at given lat lon (lat/lon is expected to be in same srs than DEM)
     * @param latitude (degree)
     * @param longitude (degree)
     * @return
     * @throws S2GeolibException
     */
    public double getElevation(double latitude, double longitude) throws RuggedException {
        double elevation = Double.NaN;
        try {
            double longitudeRad = FastMath.toRadians(longitude);
            double latitudeRad = FastMath.toRadians(latitude);

            elevation = getElevationRadian(latitudeRad, longitudeRad);
        } catch (RuggedException e) {
            throw e;
        }
        return elevation;
    }

    /**
     * Get elevation at given lat lon (lat/lon is expected to be in same SRS than DEM)
     * @param latitude (radian)
     * @param longitude (radian)
     * @return elevation (m)
     * @throws RuggedException
     */
    public double getElevationRadian(double latitudeInRad, double longitudeInRad) throws RuggedException {
        double elevation = Double.NaN;
        try {
            final double latitudeNorm = MathUtils.normalizeAngle(latitudeInRad, FastMath.PI / 2);
            final double longitudeNorm = MathUtils.normalizeAngle(longitudeInRad, 0);
            SimpleTile tile = tilesCache.getTile(latitudeNorm, longitudeNorm);
            elevation = tile.interpolateElevation(latitudeNorm, longitudeNorm);
        } catch (RuggedException e) {
            throw e;
        }
        return elevation;
    }
}









