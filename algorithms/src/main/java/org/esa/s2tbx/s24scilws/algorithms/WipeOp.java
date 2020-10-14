package org.esa.s2tbx.s24scilws.algorithms;

import org.esa.snap.core.datamodel.Band;
import org.esa.snap.core.datamodel.Mask;
import org.esa.snap.core.datamodel.Product;
import org.esa.snap.core.datamodel.ProductData;
import org.esa.snap.core.datamodel.SampleCoding;
import org.esa.snap.core.datamodel.TiePointGrid;
import org.esa.snap.core.datamodel.VirtualBand;
import org.esa.snap.core.gpf.OperatorException;
import org.esa.snap.core.gpf.OperatorSpi;
import org.esa.snap.core.gpf.annotations.OperatorMetadata;
import org.esa.snap.core.gpf.annotations.Parameter;
import org.esa.snap.core.gpf.annotations.SourceProduct;
import org.esa.snap.core.gpf.annotations.TargetProduct;
import org.esa.snap.core.gpf.pointop.PixelOperator;
import org.esa.snap.core.gpf.pointop.ProductConfigurer;
import org.esa.snap.core.gpf.pointop.Sample;
import org.esa.snap.core.gpf.pointop.SourceSampleConfigurer;
import org.esa.snap.core.gpf.pointop.TargetSampleConfigurer;
import org.esa.snap.core.gpf.pointop.WritableSample;

/**
 * The <code>MyMerisPointOp</code> serves as a template for other MERIS L1B/OLCI point ops.
 *
 * @author Norman
 */

@OperatorMetadata(
        alias = "WipeOp",
        category = "Optical/Thematic Water Processing/S2-4Sci Land and Water Study",
        version = "1.0",
        internal = false,
        description = "The 'Water Mask' operator for automatic water body detection.",
        authors = "Martino Ferrari",
        copyright = "2019")
public class WipeOp extends PixelOperator {

    private boolean hasValidPixelExpression;

    @SourceProduct(alias = "source", description = "The source product.")
    private Product sourceProduct;

    @TargetProduct
    private Product targetProduct;
    // bands 2, 3, 4, 7, 11 and 12 f
    @Parameter(label = "Blue source band",
            rasterDataNodeType = Band.class)
    private String sourceBand2;

    @Parameter(label = "Green source band",
            rasterDataNodeType = Band.class)
    // @BandParameter(minWavelength = 495, maxWavelength = 570)
    private String sourceBand3;

    @Parameter(label = "Red source band",
            rasterDataNodeType = Band.class)
    // @BandParameter(minWavelength = 600, maxWavelength = 650)
    private String sourceBand4;

    @Parameter(label = "Source band 7",
            rasterDataNodeType = Band.class)
    private String sourceBand7;

    @Parameter(label = "Source band 11",
            rasterDataNodeType = Band.class)
    private String sourceBand11;

    @Parameter(label = "Source band 12",
            rasterDataNodeType = Band.class)
    private String sourceBand12;


    /**
     * Configures all source samples that this operator requires for the computation of target samples.
     * Source sample are defined by using the provided {@link SourceSampleConfigurer}.
     * <p/>
     * <p/> The method is called by {@link #initialize()}.
     *
     * @param sampleConfigurator The configurator that defines the layout of a pixel.
     * @throws OperatorException If the source samples cannot be configured.
     */
    @Override
    protected void configureSourceSamples(SourceSampleConfigurer sampleConfigurator) throws OperatorException {
        sampleConfigurator.defineSample(0, sourceBand2);
        sampleConfigurator.defineSample(1, sourceBand3);
        sampleConfigurator.defineSample(2, sourceBand4);
        sampleConfigurator.defineSample(3, sourceBand7);
        sampleConfigurator.defineSample(4, sourceBand11);
        sampleConfigurator.defineSample(5, sourceBand12);
    }


    /**
     * Configures all target samples computed by this operator.
     * Target samples are defined by using the provided {@link TargetSampleConfigurer}.
     * <p/>
     * <p/> The method is called by {@link #initialize()}.
     *
     * @param sampleConfigurator The configurer that defines the layout of a pixel.
     * @throws OperatorException If the target samples cannot be configured.
     */
    @Override
    protected void configureTargetSamples(TargetSampleConfigurer sampleConfigurator) throws OperatorException {
        sampleConfigurator.defineSample(0, "mask");
    }

    /**
     * Configures the target product via the given {@link ProductConfigurer}. Called by {@link #initialize()}.
     * <p/>
     * Client implementations of this method usually add product components to the given target product, such as
     * {@link Band bands} to be computed by this operator,
     * {@link VirtualBand virtual bands},
     * {@link Mask masks}
     * or {@link SampleCoding sample codings}.
     * <p/>
     * The default implementation retrieves the (first) source product and copies to the target product
     * <ul>
     * <li>the start and stop time by calling {@link ProductConfigurer#copyTimeCoding()},</li>
     * <li>all tie-point grids by calling {@link ProductConfigurer#copyTiePointGrids(String...)},</li>
     * <li>the geo-coding by calling {@link ProductConfigurer#copyGeoCoding()}.</li>
     * </ul>
     * <p/>
     * Clients that require a similar behaviour in their operator shall first call the {@code super} method
     * in their implementation.
     *
     * @param productConfigurator The target product configurer.
     * @throws OperatorException If the target product cannot be configured.
     * @see Product#addBand(Band)
     * @see Product#addBand(String, String)
     * @see Product#addTiePointGrid(TiePointGrid)
     * @see Product#getMaskGroup()
     */
    @Override
    protected void configureTargetProduct(ProductConfigurer productConfigurator) {
        super.configureTargetProduct(productConfigurator);

        Product target_product = productConfigurator.getTargetProduct();
//        ProductUtils.copyBand(sourceBand2, sourceProduct, tp, true);
//        ProductUtils.copyBand(sourceBand3, sourceProduct, tp, true);
//        ProductUtils.copyBand(sourceBand4, sourceProduct, tp, true);
//        ProductUtils.copyBand(sourceBand7, sourceProduct, tp, true);
//        ProductUtils.copyBand(sourceBand11, sourceProduct, tp, true);
//        ProductUtils.copyBand(sourceBand12, sourceProduct, tp, true);
        target_product.addBand("mask", ProductData.TYPE_INT32);
    }

    /**
     * Computes the target samples from the given source samples.
     * <p/>
     * The number of source/target samples is the maximum defined sample index plus one. Source/target samples are defined
     * by using the respective sample configurator in the
     * {@link #configureSourceSamples(SourceSampleConfigurer) configureSourceSamples} and
     * {@link #configureTargetSamples(TargetSampleConfigurer) configureTargetSamples} methods.
     * Attempts to read from source samples or write to target samples at undefined sample indices will
     * cause undefined behaviour.
     *
     * @param x             The current pixel's X coordinate.
     * @param y             The current pixel's Y coordinate.
     * @param sourceSamples The source samples (= source pixel).
     * @param targetSamples The target samples (= target pixel).
     */
    @Override
    protected void computePixel(int x, int y, Sample[] sourceSamples, WritableSample[] targetSamples) {
        // Retrieve values from different bands
        double band_b = sourceSamples[0].getDouble();  //blue band
        double band_g = sourceSamples[1].getDouble();  // green band
        double band_r = sourceSamples[2].getDouble();  // red band
        double band_07 = sourceSamples[3].getDouble(); // sentinel 2 band 7
        double band_11 = sourceSamples[4].getDouble(); // sentinel 2 band 11
        double band_12 = sourceSamples[5].getDouble(); // sentinel 2 band 12

        /*
          IF:
            r OR b = 0       OR
            b11 AND b12 = 0          OR
            b11/b4 >= 0.69          OR
            b2 (blue band) < 65e-4  OR
            b11 > 0.035
          THEN mask = 0.
         */
        if ((band_b == 0 || band_r == 0) || (band_11 == 0 && band_12 == 0) || band_11 / band_r >= 0.69 || band_b < 0.0065 || band_11 > 0.035) {
            targetSamples[0].set(0);
            return;
        }

        int mask_value; // the result of the mask operation

        // Normalize rgb band values
        double norm_r = band_r / 0.15;
        double norm_g = band_g / 0.15;
        double norm_b = band_b / 0.15;

        double max_v = Math.max(Math.max(norm_r, norm_g), norm_b);  // maximum value from the rgb bands
        double min_v = Math.min(Math.min(norm_r, norm_g), norm_b);  // minimum value from the rgb bands
        double S_v = (max_v - min_v) / max_v;                       // normalization value
        if (norm_r == max_v) {
            if (max_v > -2.93 * band_11 / band_07 + 2 || band_12 / band_r > 0.28 || S_v < 0.12) {
                mask_value = 0;
            } else {
                mask_value = 1;
            }
        } else if (norm_g == max_v) {
            if (band_12 / band_r > 0.46 || S_v < 0.04) {
                mask_value = 0;
            } else if (band_r < band_b) {
                mask_value = S_v < 0.12 && max_v > 0.3 ? 0 : 1;
            } else {
                mask_value = 1;
            }
        } else {
            if (band_r - band_g > 0.001 || S_v < 0.05) {
                mask_value = 0;
            } else {
                mask_value = max_v > -1.107 * band_11 / band_b + 0.748 ? 0 : 1;
            }
        }
        targetSamples[0].set(mask_value);
    }


    public static class Spi extends OperatorSpi {
        public Spi() {
            super(WipeOp.class);
        }
    }
}

