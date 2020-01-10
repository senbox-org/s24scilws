package org.seom.dem_shadows;

/**
 * Borders min/max in latitude and max in longitudes of Earth quadrants found when creating the Geoid sub-tiles
 * @author Guylaine Prat
 */
public class EarthQuadrantsLimits {
   
   /**
    * Latitude max South (deg)
    */
   private double latMaxSouthDeg = Double.NaN;
   
   /**
    * Latitude min North (deg)
    */
   private double latMinNorthDeg = Double.NaN;

   /**
    * Longitude max FarWest (deg)
    */
   private double lonMaxFarWestDeg = Double.NaN;
   
   /**
    * Longitude max West (deg)
    */
   private double lonMaxWestDeg = Double.NaN;
   
   /**
    * Longitude max East (deg)
    */
   private double lonMaxEastDeg = Double.NaN;
   
   /**
    * Constructor
    * @param latMinNorthDeg
    * @param latMaxSouthDeg
    * @param lonMaxFarWestDeg
    * @param lonMaxWestDeg
    * @param lonMaxEastDeg
    */
   public EarthQuadrantsLimits(double latMinNorthDeg, double latMaxSouthDeg, double lonMaxFarWestDeg, double lonMaxWestDeg, double lonMaxEastDeg) {
      this.latMinNorthDeg = latMinNorthDeg;
      this.latMaxSouthDeg = latMaxSouthDeg;
      this.lonMaxFarWestDeg = lonMaxFarWestDeg;
      this.lonMaxWestDeg = lonMaxWestDeg;
      this.lonMaxEastDeg = lonMaxEastDeg;
   }

   /** Compute the latitude middle of the overlapping in North/South latitude tiles
    * in order to define the EarthQuadrant (compatible with Rugged SimpleTile.getLocation search for HAS_INTERPOLATION_NEIGHBORS)
    * @return latitude middle of the overlapping in North/South latitude tiles
    */
   public double getLimitSouthNorth(){
      return 0.5*(this.latMaxSouthDeg + this.latMinNorthDeg);
   }
   
   /**
    * Get the Longitude max FarWest (deg)
    * @return Longitude max FarWest (deg)
    */
   public double getLonMaxFarWestDeg() {
      return lonMaxFarWestDeg;
   }

   /**
    * Get the Longitude max West (deg)
    * @return Longitude max West (deg)
    */
   public double getLonMaxWestDeg() {
      return lonMaxWestDeg;
   }

   /**
    * Get the Longitude max East (deg)
    * @return Longitude max East (deg)
    */
   public double getLonMaxEastDeg() {
      return lonMaxEastDeg;
   }
}
