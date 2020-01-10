package org.seom.dem_shadows;

/**
 * Earth quadrants for Geoid tiles definition
 * @author Guylaine Prat
 *
 */
public enum EarthQuadrant {

   SOUTH_FARWEST(false, EarthHemisphere.SOUTH),
   SOUTH_WEST(false, EarthHemisphere.SOUTH),
   SOUTH_EAST(false, EarthHemisphere.SOUTH),
   SOUTH_FAREAST(true, EarthHemisphere.SOUTH),
   NORTH_FARWEST(false, EarthHemisphere.NORTH),
   NORTH_WEST(false, EarthHemisphere.NORTH),
   NORTH_EAST(false, EarthHemisphere.NORTH),
   NORTH_FAREAST(true, EarthHemisphere.NORTH);

   /**
    * Flag to tell if we are near +180 deg
    */
   private boolean flagFAREAST = false;
   
   /**
    * To tell if we are in North or South Hemisphere 
    */
   private EarthHemisphere earthHemipshereNS = null;

   
   /**
    * Constructor
    * @param isFarEast
    * @param earthHemipshereNS
    */
   private EarthQuadrant(Boolean isFarEast, EarthHemisphere earthHemipshereNS){
      this.flagFAREAST = isFarEast;
      this.earthHemipshereNS = earthHemipshereNS;
   }

   /**
    * Check if the tile is a FarEast one
    * @return
    */
   public boolean isFarEast(){
      return flagFAREAST;
   }
   
   /**
    * Get the associated hemisphere (North or South)
    * @return
    */
   public EarthHemisphere getHemisphereNS(){
      return this.earthHemipshereNS;
   }
}
