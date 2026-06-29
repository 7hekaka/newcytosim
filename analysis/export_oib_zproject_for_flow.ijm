// Export an Olympus/Bio-Formats movie as a z-projected TIFF stack for flow analysis.
//
// Argument format:
// input=/path/or/UNC/file.oib;output=/path/projected.tif;metadata=/path/meta.csv
//
// The macro opens with Bio-Formats, max-projects z slices 5-7 when available,
// falls back to the last three z slices otherwise, saves the projected time
// stack, and appends calibration/dimension metadata to a CSV.

function getArgValue(args, key) {
    parts = split(args, ";");
    prefix = key + "=";
    for (i = 0; i < parts.length; i++) {
        if (startsWith(parts[i], prefix))
            return substring(parts[i], lengthOf(prefix));
    }
    return "";
}

input = getArgValue(getArgument(), "input");
output = getArgValue(getArgument(), "output");
metadata = getArgValue(getArgument(), "metadata");

if (input == "" || output == "" || metadata == "")
    exit("Required args: input=...;output=...;metadata=...");

logFile = metadata + ".log";
File.append("start\n", logFile);
File.append("input=" + input + "\n", logFile);
File.append("output=" + output + "\n", logFile);
File.append("metadata=" + metadata + "\n", logFile);

setBatchMode(true);
File.append("before Bio-Formats Importer\n", logFile);
run("Bio-Formats Importer", "open=[" + input + "] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
File.append("after Bio-Formats Importer\n", logFile);
File.append("nImages=" + nImages + "\n", logFile);
if (nImages == 0)
    exit("Bio-Formats did not open an image for " + input);

title = getTitle();
getDimensions(width, height, channels, slices, frames);
File.append("dimensions=" + width + "," + height + "," + channels + "," + slices + "," + frames + "\n", logFile);
getVoxelSize(pixelWidth, pixelHeight, voxelDepth, unit);
File.append("voxel=" + pixelWidth + "," + pixelHeight + "," + voxelDepth + "," + unit + "\n", logFile);

requestedStart = 5;
requestedStop = 7;
if (slices >= requestedStop) {
    startSlice = requestedStart;
    stopSlice = requestedStop;
} else {
    stopSlice = slices;
    startSlice = slices - 2;
    if (startSlice < 1)
        startSlice = 1;
}

run("Z Project...", "start=" + startSlice + " stop=" + stopSlice + " projection=[Max Intensity] all");
File.append("after Z Project\n", logFile);
projectedTitle = getTitle();
saveAs("Tiff", output);
File.append("after saveAs\n", logFile);

if (!File.exists(metadata)) {
    File.append(
        "input,output,title,width,height,channels,slices,frames,pixel_width,pixel_height,voxel_depth,unit,z_start,z_stop\n",
        metadata
    );
}

File.append(
    "\"" + input + "\",\"" + output + "\",\"" + title + "\"," +
    width + "," + height + "," + channels + "," + slices + "," + frames + "," +
    pixelWidth + "," + pixelHeight + "," + voxelDepth + ",\"" + unit + "\"," +
    startSlice + "," + stopSlice + "\n",
    metadata
);

selectWindow(projectedTitle);
close();
if (isOpen(title)) {
    selectWindow(title);
    close();
}
setBatchMode(false);
File.append("done\n", logFile);
