import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.process.FloatProcessor;
import loci.formats.FormatTools;
import loci.formats.ImageReader;
import loci.formats.meta.IMetadata;
import loci.formats.MetadataTools;
import ome.units.UNITS;
import ome.units.quantity.Length;
import ome.units.quantity.Time;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class ExportOibZProject {
    private static double pixelValue(byte[] buf, int offset, int pixelType, boolean little) {
        switch (pixelType) {
            case FormatTools.UINT8:
                return buf[offset] & 0xff;
            case FormatTools.INT8:
                return buf[offset];
            case FormatTools.UINT16: {
                int b0 = buf[offset] & 0xff;
                int b1 = buf[offset + 1] & 0xff;
                return little ? (b0 | (b1 << 8)) : ((b0 << 8) | b1);
            }
            case FormatTools.INT16: {
                int b0 = buf[offset] & 0xff;
                int b1 = buf[offset + 1] & 0xff;
                short value = (short) (little ? (b0 | (b1 << 8)) : ((b0 << 8) | b1));
                return value;
            }
            case FormatTools.UINT32:
            case FormatTools.INT32: {
                long b0 = buf[offset] & 0xffL;
                long b1 = buf[offset + 1] & 0xffL;
                long b2 = buf[offset + 2] & 0xffL;
                long b3 = buf[offset + 3] & 0xffL;
                long value = little ? (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24))
                                    : ((b0 << 24) | (b1 << 16) | (b2 << 8) | b3);
                if (pixelType == FormatTools.INT32)
                    return (int) value;
                return value;
            }
            case FormatTools.FLOAT: {
                int b0 = buf[offset] & 0xff;
                int b1 = buf[offset + 1] & 0xff;
                int b2 = buf[offset + 2] & 0xff;
                int b3 = buf[offset + 3] & 0xff;
                int bits = little ? (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24))
                                  : ((b0 << 24) | (b1 << 16) | (b2 << 8) | b3);
                return Float.intBitsToFloat(bits);
            }
            default:
                throw new IllegalArgumentException("Unsupported pixel type: " + FormatTools.getPixelTypeString(pixelType));
        }
    }

    private static double lengthUm(Length value) {
        if (value == null)
            return Double.NaN;
        try {
            return value.value(UNITS.MICROMETER).doubleValue();
        } catch (Exception exc) {
            return Double.NaN;
        }
    }

    private static double timeSeconds(Time value) {
        if (value == null)
            return Double.NaN;
        try {
            return value.value(UNITS.SECOND).doubleValue();
        } catch (Exception exc) {
            return Double.NaN;
        }
    }

    private static void appendMetadata(
        String metadataCsv,
        String input,
        String output,
        int sizeX,
        int sizeY,
        int sizeZ,
        int sizeC,
        int sizeT,
        int zStartOneBased,
        int zStopOneBased,
        int pixelType,
        double pixelWidthUm,
        double pixelHeightUm,
        double voxelDepthUm,
        double timeIncrementS
    ) throws IOException {
        File file = new File(metadataCsv);
        boolean writeHeader = !file.exists();
        try (PrintWriter writer = new PrintWriter(new FileWriter(file, true))) {
            if (writeHeader) {
                writer.println("input,output,width,height,z_slices,channels,frames,z_start,z_stop,pixel_type,pixel_width_um,pixel_height_um,voxel_depth_um,time_increment_s");
            }
            writer.printf(
                "\"%s\",\"%s\",%d,%d,%d,%d,%d,%d,%d,\"%s\",%.10g,%.10g,%.10g,%.10g%n",
                input.replace("\"", "\"\""),
                output.replace("\"", "\"\""),
                sizeX,
                sizeY,
                sizeZ,
                sizeC,
                sizeT,
                zStartOneBased,
                zStopOneBased,
                FormatTools.getPixelTypeString(pixelType),
                pixelWidthUm,
                pixelHeightUm,
                voxelDepthUm,
                timeIncrementS
            );
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.err.println("Usage: ExportOibZProject <input.oib> <output.tif> <metadata.csv>");
            System.exit(2);
        }

        String input = args[0];
        String output = args[1];
        String metadataCsv = args[2];

        IMetadata metadata = MetadataTools.createOMEXMLMetadata();
        ImageReader reader = new ImageReader();
        reader.setMetadataStore(metadata);
        reader.setId(input);
        reader.setSeries(0);

        int sizeX = reader.getSizeX();
        int sizeY = reader.getSizeY();
        int sizeZ = reader.getSizeZ();
        int sizeC = reader.getSizeC();
        int sizeT = reader.getSizeT();
        int pixelType = reader.getPixelType();
        int bytesPerPixel = FormatTools.getBytesPerPixel(pixelType);
        boolean little = reader.isLittleEndian();

        int requestedStart = 5;
        int requestedStop = 7;
        int zStart = sizeZ >= requestedStop ? requestedStart - 1 : Math.max(0, sizeZ - 3);
        int zStopExclusive = sizeZ >= requestedStop ? requestedStop : sizeZ;
        if (zStart >= zStopExclusive) {
            throw new IllegalArgumentException("Invalid z projection range for sizeZ=" + sizeZ);
        }

        ImageStack stack = new ImageStack(sizeX, sizeY);
        int pixels = sizeX * sizeY;
        for (int t = 0; t < sizeT; t++) {
            float[] maxProjection = new float[pixels];
            for (int i = 0; i < pixels; i++)
                maxProjection[i] = Float.NEGATIVE_INFINITY;

            for (int z = zStart; z < zStopExclusive; z++) {
                int planeIndex = reader.getIndex(z, 0, t);
                byte[] plane = reader.openBytes(planeIndex);
                for (int i = 0; i < pixels; i++) {
                    double value = pixelValue(plane, i * bytesPerPixel, pixelType, little);
                    if (value > maxProjection[i])
                        maxProjection[i] = (float) value;
                }
            }
            stack.addSlice("t=" + (t + 1), new FloatProcessor(sizeX, sizeY, maxProjection));
        }

        ImagePlus projected = new ImagePlus(new File(input).getName() + "_zproject", stack);
        projected.setDimensions(1, 1, sizeT);
        if (sizeT > 1)
            projected.setOpenAsHyperStack(true);

        double pixelWidthUm = lengthUm(metadata.getPixelsPhysicalSizeX(0));
        double pixelHeightUm = lengthUm(metadata.getPixelsPhysicalSizeY(0));
        double voxelDepthUm = lengthUm(metadata.getPixelsPhysicalSizeZ(0));
        double timeIncrementS = timeSeconds(metadata.getPixelsTimeIncrement(0));

        Calibration cal = projected.getCalibration();
        if (!Double.isNaN(pixelWidthUm) && pixelWidthUm > 0)
            cal.pixelWidth = pixelWidthUm;
        if (!Double.isNaN(pixelHeightUm) && pixelHeightUm > 0)
            cal.pixelHeight = pixelHeightUm;
        if (!Double.isNaN(voxelDepthUm) && voxelDepthUm > 0)
            cal.pixelDepth = voxelDepthUm;
        cal.setUnit("micron");
        if (!Double.isNaN(timeIncrementS) && timeIncrementS > 0) {
            cal.frameInterval = timeIncrementS;
            cal.setTimeUnit("sec");
        }
        projected.setCalibration(cal);

        File outputFile = new File(output);
        File parent = outputFile.getParentFile();
        if (parent != null)
            parent.mkdirs();
        boolean saved = sizeT > 1
            ? new FileSaver(projected).saveAsTiffStack(output)
            : new FileSaver(projected).saveAsTiff(output);
        if (!saved)
            throw new IOException("Failed to save " + output);

        appendMetadata(
            metadataCsv,
            input,
            output,
            sizeX,
            sizeY,
            sizeZ,
            sizeC,
            sizeT,
            zStart + 1,
            zStopExclusive,
            pixelType,
            pixelWidthUm,
            pixelHeightUm,
            voxelDepthUm,
            timeIncrementS
        );

        reader.close();
    }
}
