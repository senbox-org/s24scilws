package org.seom.SunIntersectionGrid.createGrid;

import java.io.File;
import java.io.IOException;
import java.util.Locale;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import org.gdal.gdal.Dataset;
import org.gdal.gdal.Driver;
import org.gdal.gdal.gdal;
import org.gdal.gdalconst.gdalconst;
import org.gdal.osr.CoordinateTransformation;
import org.gdal.osr.SpatialReference;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.data.DataProvidersManager;
import org.orekit.data.DirectoryCrawler;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.rugged.errors.RuggedException;
import org.orekit.rugged.intersection.IntersectionAlgorithm;
import org.orekit.rugged.intersection.duvenhage.DuvenhageAlgorithm;
import org.orekit.rugged.utils.ExtendedEllipsoid;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.TimeStampedPVCoordinates;
import org.seom.dem_shadows.DemFileManager;
import org.seom.dem_shadows.DemManager;
import org.seom.dem_shadows.GeoidManager;
import org.seom.dem_shadows.SRTM1sFileManager;
import org.seom.utilities.ComputeUtils;
//import org.seom.SunIntersectionGrid.src.fileReader.ElementReader;
import org.seom.utilities.RuggedConfiguration;



public class createTopoShadowGrid {


    public static void main(String[] args) {

        if(args.length < 5)
        {
            System.out.println("usage : RuggedConfigFile RefTime OutputGridFilename Grid_ULX Grid_ULY Grid_StepX Grid_StepY Grid_SizeX GridSizeY (EPSG)");
            System.out.println("By default EPSG is 4326 (WGS84)");
            return;
        }
        System.out.println("Progress[%]: 0.0");

        String confFile = args[0];
        RuggedConfiguration config;
        try {
            config = new RuggedConfiguration(confFile);
            File orekitData = new File(config.getString(RuggedConfiguration.OREKITDATA));
            DataProvidersManager.getInstance().addProvider(new DirectoryCrawler(orekitData));


            String demRootDir = config.getString(RuggedConfiguration.DEMROOTDIR);
            String geoidPath = config.getString(RuggedConfiguration.GEOIDPATH);

            TimeScale utc = TimeScalesFactory.getUTC();
            AbsoluteDate referenceDate = new AbsoluteDate(args[1], utc);
            // Reference frames
            Frame itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, false);

            // Sun's definition
            CelestialBody sun = CelestialBodyFactory.getSun();

            // Create Earth body
            ExtendedEllipsoid earth = new ExtendedEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                                            Constants.WGS84_EARTH_FLATTENING,
                                                            itrf);

            // Initialize Geoid manager
            GeoidManager geoidManager = GeoidManager.initGeoidManager(geoidPath);

            // Initialize the DEM manager
            System.out.println("Initializing DEM manager...");
            DemFileManager demFileManager = new SRTM1sFileManager(demRootDir);
            System.out.format("SRTM root dir %s %n", demFileManager.getDemRootDir());

            // tests of correct file are performed in DemFileManager inherited classes
            demFileManager.findRasterFile();
            DemManager demManager = DemManager.initDemManager(demFileManager);

            // Duvenhage algorithm implementation
            IntersectionAlgorithm duvenhage = new DuvenhageAlgorithm(demManager, 4, false);

            // Get Sun's position in Earth frame (ITRF) for image reference date
            TimeStampedPVCoordinates sunPV = sun.getPVCoordinates(referenceDate, itrf);
            Vector3D sunPosition = sunPV.getPosition();
            System.out.format(Locale.US, "sunPosition: X=%8.3f, Y = %8.3f, Z = %8.3f%n", sunPosition.getX(), sunPosition.getY(), sunPosition.getZ());

            //final double sunAzimuth = localFrame.getAzimuth(sunPosition, itrf, referenceDate);
            //final double sunElevation = localFrame.getElevation(sunPosition, itrf, referenceDate);
            //System.out.format(Locale.US, "sunAzimuth  = %3.1f ° sunZenith = %3.1f °",FastMath.toDegrees(sunAzimuth), 90.0 - FastMath.toDegrees(sunElevation));


            //loop over the grid
            double xUL = new Double(args[3]);
            double yUL = new Double(args[4]);
            double xRes  = new Double(args[5]);
            double yRes  = new Double(args[6]);
            int xSize = new Integer(args[7]);
            int ySize = new Integer(args[8]);

            int EPSG = 4326;
            if(args.length == 10)
            {
                EPSG=new Integer(args[9]);
            }
            int gridSize = xSize * ySize;
            // Create a grid with respect of ULLR and spacing path
            double[] xGrid = ComputeUtils.gridCreation(xUL, xRes, xSize);
            double[] yGrid = ComputeUtils.gridCreation(yUL, yRes, ySize);
            float[] outGrid = new float[gridSize];


            gdal.AllRegister();

            Driver driver = gdal.GetDriverByName("GTiff");
            Dataset dataset = driver.Create(args[2], xSize, ySize, 1, gdalconst.GDT_Float32);
            // Use the sensor grid to define both grids (on board and on ground ...)
            CoordinateTransformation ct =  null;
            if(EPSG != 4326)
            {
                SpatialReference src = new SpatialReference("");
                src.ImportFromEPSG(EPSG);
                SpatialReference dst = new SpatialReference("");
                dst.ImportFromEPSG(4326);
                ct = new CoordinateTransformation(src, dst);
            }
            DecimalFormat formatter = new DecimalFormat("#0.0", new DecimalFormatSymbols(Locale.US));
            int currentPurcentage = 0;
            for (int line=0; line < ySize; line++) {
                for (int col=0; col < xSize; col++) {
                    double x = xGrid[col];
                    double y = yGrid[line];
                    //System.out.println("Pixel_deg : "+pixel_deg+" Line_deg : "+line_deg);
                    // Conversion EPSG_origin to EPSG:4326
                    double altitudeA = 0.0;

                    double[] point = {x, y,altitudeA};
                    if (ct != null)
                        ct.TransformPoint(point);

                    double latitudeA = FastMath.toRadians(point[1]);
                    double longitudeA = FastMath.toRadians(point[0]);
                    altitudeA = duvenhage.getElevation(latitudeA, longitudeA);

                    GeodeticPoint gpA = new GeodeticPoint(latitudeA, longitudeA, altitudeA);
                    //get sun zenith/azimuth wrt to point A

                    Vector3D pA = earth.transform(gpA);
                    //            System.out.format(Locale.US, "pA: X=%8.3f, Y = %8.3f, Z = %8.3f%n", pA.getX(), pA.getY(), pA.getZ());

                    // Line of sight estimation
                    Vector3D los = pA.subtract(sunPosition).normalize();

                    // Get B geodetic point position
                    GeodeticPoint gpB = duvenhage.intersection(earth, sunPosition, los);
                    //GeodeticPoint gpB =  gpA;
                    //System.out.format(Locale.US, "Ground point B: φ = %8.3f °, λ = %8.3f °, h = %8.3f m%n",
                    //                  FastMath.toDegrees(gpB.getLatitude()), FastMath.toDegrees(gpB.getLongitude()), gpB.getAltitude());

                    Vector3D pB = earth.transform(gpB);
                    //            System.out.format(Locale.US, "pB: X=%8.3f, Y = %8.3f, Z = %8.3f%n", pB.getX(), pB.getY(), pB.getZ());



                    outGrid[col+line*xSize] = (float) Vector3D.distance(pA, pB);
                }
                int comparePurcentage=(int)(100*line/(float)(ySize+2));
                System.out.println(comparePurcentage);
                if(currentPurcentage<comparePurcentage){
                    currentPurcentage = comparePurcentage;
                    System.out.println("Progress[%]: "+formatter.format(currentPurcentage));
                }
            }

            double [] geotransform = {xUL, xRes, 0, yUL, 0, yRes};
            dataset.SetGeoTransform(geotransform);
            SpatialReference srs = new SpatialReference();
            srs.ImportFromEPSG(EPSG);
            dataset.SetProjection(srs.ExportToPrettyWkt());
            dataset.GetRasterBand(1).WriteRaster(0, 0, xSize, ySize, xSize, ySize, gdalconst.GDT_Float32, outGrid);


            dataset.FlushCache();
            dataset.delete();
            System.out.println("Progress[%]: 100");
            System.out.println("Created file: "+args[2]);

        }  catch (IOException e) {
            e.printStackTrace();
        } catch (OrekitException e) {
            e.printStackTrace();
        } catch (RuggedException e) {
            e.printStackTrace();
        }
    }
}


