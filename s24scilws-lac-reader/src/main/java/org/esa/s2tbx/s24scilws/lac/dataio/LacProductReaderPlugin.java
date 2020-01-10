package org.esa.s2tbx.s24scilws.lac.dataio;

import org.esa.s2tbx.dataio.readers.BaseProductReaderPlugIn;
import org.esa.s2tbx.dataio.readers.PathUtils;
import org.esa.snap.core.dataio.DecodeQualification;
import org.esa.snap.core.dataio.ProductReader;
import org.esa.snap.core.dataio.ProductReaderPlugIn;
import org.esa.snap.core.datamodel.RGBImageProfile;
import org.esa.snap.core.datamodel.RGBImageProfileManager;
import org.esa.snap.core.util.io.SnapFileFilter;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

/**
 * Created by obarrile on 25/12/2019.
 */
public class LacProductReaderPlugin implements ProductReaderPlugIn {

    @Override
    public DecodeQualification getDecodeQualification(Object input) {
        Path imageIOInputPath = convertInputToPath(input);
        if(imageIOInputPath.endsWith("output_DATA")) {
            return DecodeQualification.INTENDED;
        }
        if(imageIOInputPath.endsWith("MeanO3.txt")) {
            return DecodeQualification.INTENDED;
        }
        return DecodeQualification.UNABLE;
    }

    @Override
    public Class[] getInputTypes() {
        return LacReaderConstants.LAC_READER_INPUT_TYPES;
    }

    @Override
    public ProductReader createReaderInstance() {
        return new LacProductReader(this);
    }

    @Override
    public String[] getFormatNames() {
        return LacReaderConstants.LAC_FORMAT_NAMES;
    }

    @Override
    public String[] getDefaultFileExtensions() {
        return LacReaderConstants.LAC_DEFAULT_EXTENSIONS;
    }

    @Override
    public String getDescription(Locale locale) {
        return LacReaderConstants.LAC_PRODUCT_DESCRIPTION;
    }

    @Override
    public SnapFileFilter getProductFileFilter() {
        return new SnapFileFilter(getFormatNames()[0], getDefaultFileExtensions(), getDescription(null));
    }

    public static Path convertInputToPath(Object input) {
        if (input == null) {
            throw new NullPointerException();
        } else if (input instanceof File) {
            return ((File) input).toPath();
        } else if (input instanceof Path) {
            return (Path) input;
        } else if (input instanceof String) {
            return Paths.get((String) input);
        } else {
            throw new IllegalArgumentException("Unknown input '" + input + "'.");
        }
    }

}
