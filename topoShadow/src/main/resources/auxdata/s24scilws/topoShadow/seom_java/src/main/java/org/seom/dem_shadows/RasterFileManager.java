package org.seom.dem_shadows;
import java.io.File;

import org.hipparchus.exception.DummyLocalizable;
import org.orekit.rugged.errors.RuggedException;

public class RasterFileManager extends DemFileManager {

    protected String demfileName;

    /**
     * {@inheritDoc}
     */
    public RasterFileManager(String demFilename) {
        super(demFilename);
        demfileName = demFilename;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected boolean findRasterFile(String fileName) throws RuggedException {
        File demFile = new File(fileName);

        if (!demFile.exists()) {
            String str = "raster not found : " + fileName;
            DummyLocalizable message = new DummyLocalizable(str);
            throw new RuggedException(message);
        }
        demfileName = fileName;
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getRasterFilePath(double latitude, double longitude) {
        return demfileName;
    }

}
