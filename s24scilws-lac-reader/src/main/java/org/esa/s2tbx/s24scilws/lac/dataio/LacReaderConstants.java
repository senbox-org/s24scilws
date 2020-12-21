package org.esa.s2tbx.s24scilws.lac.dataio;

import java.io.File;

/**
 * Created by obarrile on 25/12/2019.
 */
public class LacReaderConstants {

    // constants for plugins
    public static final Class[] LAC_READER_INPUT_TYPES = new Class[]{String.class, File.class};
    public static final String LAC_PRODUCT_DESCRIPTION = "LAC processor Data Products";
    public static final String[] LAC_DEFAULT_EXTENSIONS = new String[]{""};
    public static final String[] LAC_FORMAT_NAMES = new String[]{"LAC"};
    public static final String[] MINIMAL_PRODUCT_PATTERNS = new String[] {
            "MTD_LAC.xml"
    };
}
