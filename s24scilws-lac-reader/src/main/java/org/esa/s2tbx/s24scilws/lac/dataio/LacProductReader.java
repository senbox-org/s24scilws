package org.esa.s2tbx.s24scilws.lac.dataio;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.s2tbx.dataio.VirtualDirEx;
import org.esa.s2tbx.dataio.readers.BaseProductReaderPlugIn;
import org.esa.s2tbx.dataio.s2.S2BandConstants;
import org.esa.snap.core.dataio.AbstractProductReader;
import org.esa.snap.core.dataio.ProductReaderPlugIn;
import org.esa.snap.core.datamodel.Band;
import org.esa.snap.core.datamodel.Product;
import org.esa.snap.core.datamodel.ProductData;
import org.esa.snap.core.util.ProductUtils;
import org.esa.snap.dataio.geotiff.GeoTiffProductReader;
import org.esa.snap.dataio.geotiff.GeoTiffProductReaderPlugIn;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Arrays;  

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.DocumentBuilder;
import org.w3c.dom.Document;
import org.xml.sax.SAXException;
import org.w3c.dom.Node;

/**
 * Created by obarrile on 25/12/2019.
 */
public class LacProductReader extends AbstractProductReader {

    private static final Logger logger = Logger.getLogger(LacProductReader.class.getName());
    private VirtualDirEx virtualDir;
    private List<Product> associatedProducts;

    protected LacProductReader(ProductReaderPlugIn readerPlugIn) {
        super(readerPlugIn);
        this.associatedProducts = new ArrayList<>();
    }

    private static S2BandConstants getBandFromName(String name) {
        for (S2BandConstants band : S2BandConstants.values()) {
            if (band.getFilenameBandId().equals(name))
                return band;
        }
        return null;
    }

    @Override
    protected Product readProductNodesImpl() throws IOException {

        Path inputPath = BaseProductReaderPlugIn.convertInputToPath(super.getInput());
        Path metadaPath = inputPath.normalize();
        if (inputPath.getFileName().toString().equals("MTD_LAC.xml")) {
            inputPath.resolveSibling("GRANULE");
        }

        if (logger.isLoggable(Level.FINE)) {
            logger.log(Level.FINE, "Reading LAC product from the file '" + inputPath.toString() + "'.");
        }

        this.virtualDir = VirtualDirEx.build(inputPath, false, true);
        if (this.virtualDir == null) {
            throw new NullPointerException("The virtual dir is null for input path '" + inputPath.toString() + "'.");
        }

        // read all Tifs
        ArrayList<Band> bandList = new ArrayList<>();
        Band referenceBand = null; // THis will be the largest band
        String[] files = this.virtualDir.listAll();
        Arrays.sort(files);
        for (String file : files) {
            if (file.endsWith(".TIF")) {
                Band band = readGeoTiffProductBand(file, 0); // all TIFs in this product should have only one band
                bandList.add(band);
                if (referenceBand == null || referenceBand.getRasterWidth() < band.getRasterWidth()) {
                    referenceBand = band;
                }
            }
        }
        File file = metadaPath.toFile();
        // an instance of factory that gives a document builder
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        // an instance of builder to parse the specified xml file
        DocumentBuilder db;
        try {
            db = dbf.newDocumentBuilder();
            Document doc = db.parse(file);  
            doc.getDocumentElement().normalize();  
            Node node = doc.getElementsByTagName("PRODUCT_URI").item(0);  
            //create Product
            Product product = new Product(node.getTextContent(), "LAC", referenceBand.getRasterWidth(), referenceBand.getRasterHeight());
            product.setDescription("LAC output");
            // set File Location
            product.setFileLocation(inputPath.toFile());

            product.setSceneGeoCoding(referenceBand.getGeoCoding());
            product.setAutoGrouping("Bands");

            //addBands
            for(Band srcBand : bandList) {
                String file_name = srcBand.getProduct().getName();
                String band_name = file_name;
                boolean visual_band = false;
                S2BandConstants bandInfo = null;
                if (file_name.startsWith("O3-")){
                    band_name = "LAC_O3";
                } else if (file_name.startsWith("C_")){
                    band_name = file_name.split("_")[3];
                    bandInfo = LacProductReader.getBandFromName(band_name);
                    band_name = bandInfo.getPhysicalName();
                    visual_band = true;

                } else {
                    band_name = "LAC_"+file_name.split("_")[0];
                }
                Band targetBand = new Band(band_name, srcBand.getDataType(), srcBand.getRasterWidth(), srcBand.getRasterHeight());
                product.addBand(targetBand);
                ProductUtils.copyGeoCoding(srcBand, targetBand);
                targetBand.setNoDataValue(srcBand.getNoDataValue());
                targetBand.setNoDataValueUsed(true);
                if (visual_band) {
                    targetBand.setSpectralWavelength((float)bandInfo.getWavelengthCentral());
                    targetBand.setSpectralBandwidth((float)(bandInfo.getWavelengthMax() - bandInfo.getWavelengthMin()));
                }
                //targetBand.setScalingFactor(1.0d);
                //targetBand.setScalingOffset(0.0d);
                targetBand.setSampleCoding(srcBand.getSampleCoding());
                targetBand.setImageInfo(srcBand.getImageInfo());
                //targetBand.setDescription(""); //TODO
                targetBand.setSourceImage(srcBand.getSourceImage());
            }

            return product;
        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            e.printStackTrace();
        }
        return null;
    }

    @Override
    protected void readBandRasterDataImpl(int sourceOffsetX, int sourceOffsetY, int sourceWidth, int sourceHeight, int sourceStepX, int sourceStepY, Band destBand, int destOffsetX, int destOffsetY, int destWidth, int destHeight, ProductData destBuffer, ProgressMonitor pm) throws IOException {
        // do nothing
    }


    // get the bands and include the product in associated product, to be properly closed when closing Muscate product
    private Band readGeoTiffProductBand(String pathString, int bandIndex) {
        Band band = null;
        try {
            File inputFile;
            try {
                inputFile = this.virtualDir.getFile(pathString);
            } catch (FileNotFoundException e) {
                String fileName = pathString.substring(pathString.lastIndexOf("/") + 1);
                inputFile = this.virtualDir.getFile(fileName);
            }

            GeoTiffProductReaderPlugIn geoTiffReaderPlugIn = new GeoTiffProductReaderPlugIn();
            GeoTiffProductReader geoTiffProductReader = new GeoTiffProductReader(geoTiffReaderPlugIn);
            Product tiffProduct = geoTiffProductReader.readProductNodes(inputFile, null);
            this.associatedProducts.add(tiffProduct);
            band = tiffProduct.getBandAt(bandIndex);
        } catch (IOException e) {
            logger.warning(String.format("Unable to get band %d of the product: %s", bandIndex, pathString));
        }
        return band;
    }

    @Override
    public void close() throws IOException {
        super.close();

        for (Product product : this.associatedProducts) {
            product.dispose();
        }
        this.associatedProducts.clear();
        this.virtualDir.close();
        this.associatedProducts = null;
    }
}
