package org.seom.dem_shadows;

import java.util.Locale;

import org.hipparchus.util.FastMath;

/**
 * @author gprat
 *
 */
public class Constants {

    //TODO check methods and javadoc

    // ==================================================================
    // Default values for s2geolib.properties parameters if not present
    // ==================================================================
    /**
     * Tell if rugged must enable light time correction
     * for compensating (TRUE, default) or not (FALSE) light time between ground and spacecraft.
     * Compensating this delay improves localization accuracy and is enabled by default.
     * Not compensating it is mainly useful for validation purposes against system that do not compensate it.
     **/
    public static final boolean LIGHT_TIME_CORRECTION = true;

    /**
     * Tell if rugged must enable aberration of light correction
     * for compensating (TRUE, default) or not (FALSE) aberration of light,
     * which is velocity composition between light and spacecraft when the
     * light from ground points reaches the sensor.
     * Compensating this velocity composition improves localization
     * accuracy and is enabled by default. Not compensating it is useful
     * in two cases:
     * for validation purposes against system that do not compensate it
     * or when the pixels line of sight already include the correction.
     */
    public static final boolean ABERRATION_OF_LIGHT_CORRECTION = true;


    /**
     * For velocity computation without kinematic effects (case of some CNES data)
     * By default, real data take into account kinematic effects (flag set to FALSE by default !)
     */
    public static final boolean DUMMY_PV_TRANSFORM = false;

    /**
     * When printing direct location results, tell if longitude and latitude must be written
     * if true : deg/min/sec (example: West 135° 7′ 3.33595″)
     * if false : in decimal degrees (example: -135.11759332)
     */
    public static final boolean PRINT_IN_DMS = false;

    /**
     * Default local for messages
     */
    public static final Locale DEFAULT_LOCALE = Locale.ENGLISH;


    /**
     * Earth radius in cms
     */
    public static final double EARTH_RADIUS = 637100000d;


    /**
     * One degree value in radian
     */
    public static final double ONE_DEGREE_IN_RADIAN = FastMath.toRadians(1d);

    /**
     * Default output directory for colocation grid
     */
    public static final String DEFAULT_COLOCATION_GRID_DIR_NAME = "COLOCATION_GRID";

    /**
     * Default output directory for direct loc grid
     */
    public static final String DEFAULT_DIRECT_LOCATION_GRID_DIR_NAME = "DIRECT_LOCATION_GRID";

    /**
     * Default output directory for inverse loc grid
     */
    public static final String DEFAULT_INVERSE_LOCATION_GRID_DIR_NAME = "INVERSE_LOCATION_GRID";
    /**
     * Default no data value
     */
    public static final double NO_DATA_VALUE = -99999d;



    //============================================
    // Configuration for Rugged
    //============================================

    //    /**
    //     * Full resolution step for rugged (given to RuggedBuilder.setTimeSpan)
    //     */
    //    public static final double RUGGED_FULL_RES_STEP = 0.001d;

    /**
     * Low resolution step for rugged (given to RuggedBuilder.setTimeSpan)
     */
    public static final double RUGGED_LOW_RES_STEP = 0.1d;

    /**
     * tolerance in seconds allowed for minDate and maxDate overshooting (given to RuggedBuilder.setTimeSpan)
     */
    public static final int OVERSHOOT_TOLERANCE = 10;

    /**
     * Nb cached tiles send to init rugged (given to RuggedBuilder.setDigitalElevationModel)
     */
    public static final int NB_CACHE_TILE = 8;

    /**
     * Nb time rugged ask exactly same DEM tile: means that infinite loop due to DEM tiles not overlapping; use in DEM updateTile method
     */
    public static final int NB_TIME_ASK_SAME_TILE = 1000;

    /**
     * position/velocity interpolation order to init rugged (number of points to use for position/velocity  interpolation)
     */
    public static final int PV_INTERPOLATION_ORDER = 6;

    /**
     * quaternion interpolation order to init rugged (number of points to use for attitude interpolation)
     */
    public static final int A_INTERPOLATION_ORDER = 8;


}
