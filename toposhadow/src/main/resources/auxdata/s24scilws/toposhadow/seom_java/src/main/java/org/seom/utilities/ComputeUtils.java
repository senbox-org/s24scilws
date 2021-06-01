package org.seom.utilities;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.bodies.GeodeticPoint;

/**
 * Utilities methods
 * @author CS team
 * TBN: Excerpt from Rugged tool.
 */
public class ComputeUtils {

    /** Compute an (approximate) geodetic distance in meters between geodetic points (long1, lat1) and (long2, lat2).
     * TBN: Orekit does not have such a method
     * @param earthRadius Earth radius (m)
     * @param long1 longitude of first geodetic point (rad)
     * @param lat1 latitude of first geodetic point (rad)
     * @param long2 longitude of second geodetic point (rad)
     * @param lat2 latitude of second geodetic point (rad)
     * @return distance in meters
     */
    public static double computeDistanceInMeter(double earthRadius, final double long1, final double lat1,
                                                final double long2, final double lat2) {

        // get vectors on unit sphere from angular coordinates
        final Vector3D p1 = new Vector3D(long1, lat1);
        final Vector3D p2 = new Vector3D(long2, lat2);
        return earthRadius * Vector3D.angle(p1, p2);
    }

    /** Compute an (approximate) geodetic distance in meters between two geodetic points
     * TBN: Orekit does not have such a method
     * @param earthRadius Earth radius (m)
     * @param gp1 first geodetic point
     * @param gp2 second geodetic point
     * @return distance in meters
     */
    public static double computeDistanceInMeter(double earthRadius, final GeodeticPoint gp1, final GeodeticPoint gp2) {

        return computeDistanceInMeter(earthRadius, gp1.getLongitude(), gp1.getLatitude(), gp2.getLongitude(), gp2.getLatitude());
    }


    public static double[] gridCreation(final double UL, final double res, final int n) {
        final double[] grid = new double[n];
        for (int i = 0; i < n; ++i) {
            grid[i] = UL + (i +0.5) * (res);
        }
        return grid;
    }
}
