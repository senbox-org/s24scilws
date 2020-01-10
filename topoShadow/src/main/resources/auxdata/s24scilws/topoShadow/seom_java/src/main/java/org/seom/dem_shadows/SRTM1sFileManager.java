package org.seom.dem_shadows;
import java.io.File;
import java.nio.file.DirectoryStream;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;

import org.hipparchus.exception.DummyLocalizable;
import org.hipparchus.util.FastMath;
import org.orekit.rugged.errors.RuggedException;

public class SRTM1sFileManager extends DemFileManager {

    /**
     * {@inheritDoc}
     */
    public SRTM1sFileManager(String demRootDir) {
        super(demRootDir);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean findRasterFile(String directory) throws RuggedException {
        try {
            Path dir = FileSystems.getDefault().getPath(directory);
            DirectoryStream<Path> stream = Files.newDirectoryStream(dir);
            boolean found = false;
            for (Path path : stream) {
                if (!found) {
                    File currentFile = path.toFile();
                    if (currentFile.isDirectory()) {
                        found = findRasterFile(currentFile.getAbsolutePath());
                    } else {
                        String filePath = currentFile.getAbsolutePath();
                        if (filePath.matches(".*.tif")) {
                            found = true;
                        }
                    }
                    if (found) {
                        stream.close();
                        return true;
                    }
                }
            }
            stream.close();

            String str = "raster not found  in" + demRootDir;
            DummyLocalizable message = new DummyLocalizable(str);
            throw new RuggedException(message);

        } catch (Exception e) {
            String str = "dir not found " + demRootDir;
            DummyLocalizable message = new DummyLocalizable(str);
            throw new RuggedException(message);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getRasterFilePath(double latitude, double longitude) {
        double latFloor = FastMath.floor(FastMath.toDegrees(latitude));
        double lonFloor = FastMath.floor(FastMath.toDegrees(longitude));
        int lat = (int) latFloor;
        int lon = (int) lonFloor;

        String str = "";
        String demDir = null;
        // Compute file path name :
        if (lon < 0) {
            str += "w";
        } else {
            str += "e";
        }
        str += String.format("%03d", FastMath.abs(lon));
        demDir = str;
        if (lat < 0) {
            str += "s";
        } else {
            str += "n";
        }
        String filePath = demRootDir + File.separator + demDir + File.separator + str + String.format("%02d", FastMath.abs(lat)) + ".tif";
        return filePath;
    }

}
