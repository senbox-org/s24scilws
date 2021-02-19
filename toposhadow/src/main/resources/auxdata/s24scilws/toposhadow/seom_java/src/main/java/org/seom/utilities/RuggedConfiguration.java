package org.seom.utilities;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * ADSrugged configuration file reader.
 * @author CS team
 */
public class RuggedConfiguration {

    private Properties config = null;
    private final String confFilename;

    // Configuration

    /** Value for grid sampling in degrees */
    public static final String DEMROOTDIR="dem.root.directory";
    public static final String GEOIDPATH="geoid.path";
    public static final String OREKITDATA="physical.data.dir";


    /**
     * Constructor with the conf file name
     * @param confFileName
     * @throws Exception
     */
    public RuggedConfiguration(final String confFileName) throws IOException {

        this.confFilename = confFileName;

        File confFile = new File(confFileName); // should not be null !
        InputStream fileDir = null;
        try {
            fileDir = new FileInputStream(confFile);
            config = new Properties();
            config.load(fileDir);
        } catch (IOException e) {
            System.err.println("The configuration file " + confFileName + " does not exist");
            throw e;
        }
    }

    /**
     * Get a string
     */
    public String getString(String key) {
        if (config.containsKey(key)) return config.getProperty(key);
        System.err.println("The key " + key + " in the configuration file "+ this.confFilename + " does not exist");
        return null;
    }

    /**
     * Get a string or its default value if does not exist
     */
    public String getString(String key, String default_val) {
        if (config.containsKey(key)) return config.getProperty(key);
        return default_val;
    }

    /**
     * Get a boolean
     */
    public Boolean getBoolean(String key) {
        if (config.containsKey(key)) {
            String property = config.getProperty(key);
            if (property != null && property.length() > 0) return Boolean.parseBoolean(property);
        }
        System.err.println("The key " + key + " in the configuration file "+ this.confFilename + " does not exist");
        return null;
    }
    /**
     * Get a boolean or its default value if does not exist
     */
    public Boolean getBoolean(String key, Boolean default_val) {
        if (config.containsKey(key)) {
            String property = config.getProperty(key);
            if (property != null && property.length() > 0) return Boolean.parseBoolean(property);
        }
        return default_val;
    }

    /**
     * Get a double
     */
    public Double getDouble(String key) {
        if (config.containsKey(key)) {
            String property = config.getProperty(key);
            if (property != null && property.length() > 0) {
                return Double.parseDouble(property);
            }
        }
        System.err.println("The key " + key + " in the configuration file "+ this.confFilename + " does not exist");
        return null;
    }

    /**
     * Get a double or its default value if does not exist
     */
    public Double getDouble(String key, Double default_val) {
        if (config.containsKey(key)) {
            String property = config.getProperty(key);

            if (property != null && property.length() > 0) {
                Double readValue = Double.parseDouble(property);
                if( !Double.isNaN(readValue)) {
                    return readValue;
                } else {
                    return default_val;
                }
            }
        }
        return default_val;
    }

    /**
     * Get an integer
     */
    public Integer getInteger(String key) {
        if (config.containsKey(key)) {
            String property = config.getProperty(key);

            if (property != null && property.length() > 0) {
                return Integer.parseInt(property);
            }
        }
        System.err.println("The key " + key + " in the configuration file "+ this.confFilename + " does not exist");
        return null;
    }

    /**
     * Get an integer or its default value if does not exist
     */
    public Integer getInteger(String key, Integer default_val) {
        if (config.containsKey(key)) {
            String property = config.getProperty(key);
            if (property != null && property.length() > 0) {
                return Integer.parseInt(property);
            }
        }
        return default_val;
    }
}
