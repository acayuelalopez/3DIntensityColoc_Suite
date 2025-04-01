import ij.IJ
import ij.ImagePlus
import ij.ImageStack
import ij.WindowManager
import ij.measure.Calibration
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import inra.ijpb.label.LabelImages
import inra.ijpb.measure.region3d.RegionAnalyzer3D
import inra.ijpb.morphology.Strel
import loci.plugins.BF
import loci.plugins.in.ImporterOptions
import mcib3d.geom.Objects3DPopulation
import mcib3d.image3d.ImageInt
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest
import java.util.stream.Stream;

// INPUT UI
//
#
@File(label = "Input File Directory", style = "directory") inputFilesRef
#
@File(label = "Output directory", style = "directory") outputDir
#
@File(label = "Green model", style = "file") greenModel
#
@File(label = "Red model", style = "file") redModel
#
@Integer(label = "Reference Channel", value = 1) refIndex
#
@Integer(label = "Target Channel", value = 2) targetIndex


// IDE
//
//
//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -inputFileDir Ref: " + inputFilesRef)
IJ.log("    -outputDir: " + outputDir)
IJ.log("                                                           ");
/** Get files (images) from input directory */
def listOfFilesRef = inputFilesRef.listFiles(); ;
/** Create table for KST pvalue */
def tablePerFile = new ResultsTable()
/** Create table for overlap measure (all labels,per image) */
def totalTable = new ResultsTable();
/** Define wt/ko clones */
def clones = new String[]{"WT21", "KO8", "WT28", "KO28", "WT6", "KO19", "WT10", "KO2", "WT18", "KO20", "WT30", "KO25"}

for (def i = 0; i < listOfFilesRef.length; i++) {

    /** Create table for overlap measure (per label) */
    //def rb = new ResultsBuilder();
    def tableOverlap = new ResultsTable();

    //if (!listOfFilesRef[i].getName().contains("DS") & !listOfFilesRef[i].getName().contains("_seg.npy")) {

//        /** Importer options for .lif file */
//        def options = new ImporterOptions();
//        options.setId(inputFilesDir.getAbsolutePath() + File.separator + listOfFiles[i].getName());
//        options.setSplitChannels(false);
//        options.setSplitTimepoints(false);
//        options.setSplitFocalPlanes(false);
//        options.setAutoscale(true);
//        options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
//        options.setStackOrder(ImporterOptions.ORDER_XYCZT);
//        options.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE);
//        options.setCrop(false);
//        options.setOpenAllSeries(true);
//        def imps = BF.openImagePlus(options);
//        IJ.log("Analyzing file: " + listOfFiles[i].getName());


    //for (def j = 0.intValue(); j < imps.length; j++) {
    //IJ.log("   -Analyzing serie: " + (j + 1));
    /** Get image serie per lif */
//            def imp = imps[j]
//            /** Get channels separately */
//            def channels = ChannelSplitter.split(imp)
//            def blueCh = channels[0]
//            def refCh = channels[refIndex.intValue()]
//            def targetCh = channels[targetIndex.intValue()]
//        def refCh = new ImagePlus(listOfFilesRef[i].getAbsolutePath())
//        def targetCh = new ImagePlus(listOfFilesTarget[i].getAbsolutePath())
    def refChSeg = new ImagePlus(listOfFilesRef[i].getAbsolutePath())
    def targetChSeg = new ImagePlus(listOfFilesRef[i].getAbsolutePath().replaceAll("green", "red").replaceAll("G", "R").replaceAll("RRO", "RGO"))


//        if (refIndex == 1) {


//                /** Analyzing ref channel (signal) == green */
//                refCh.show()
//                IJ.run(refCh, "Cellpose Advanced (custom model)", String.format("diameter=0 cellproba_threshold=0.0 flow_threshold=0.4 anisotropy=1.0 diam_threshold=9999999 model_path=%s model=%s nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=", greenModel.getAbsolutePath(), greenModel.getAbsolutePath()));
//                refCh.hide()
//                refChSeg = WindowManager.getImage(refCh.getShortTitle() + "-cellpose")
//                IJ.saveAs(refChSeg, "Tiff", outputDir.getAbsolutePath() + File.separator + imps[j].getTitle().replaceAll("/", ""))
//                refChSeg.hide()
//
//                /** Analyzing target channel (signal) == red */
//                targetCh.show()
//                IJ.run(targetCh, "Cellpose Advanced (custom model)", String.format("diameter=0 cellproba_threshold=0.0 flow_threshold=1.5 anisotropy=1.0 diam_threshold=9999999 model_path=%s model=%s nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=", redModel.getAbsolutePath(), redModel.getAbsolutePath()));
//                targetCh.hide()
//                targetChSeg = WindowManager.getImage(targetCh.getShortTitle() + "-cellpose")
//                IJ.saveAs(targetChSeg, "Tiff", outputDir.getAbsolutePath() + File.separator + imps[j].getTitle().replaceAll("/", ""))
//                targetChSeg.hide()

//        } else {
//                /** Analyzing ref channel (signal) == red */
//                IJ.run(refCh, "Cellpose Advanced (custom model)", String.format("diameter=0 cellproba_threshold=0.0 flow_threshold=0.4 anisotropy=1.0 diam_threshold=9999999 model_path=%s model=%s nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=", redModel.getAbsolutePath(), redModel.getAbsolutePath()));
//                refCh.hide()
//                refChSeg = WindowManager.getImage(refCh.getShortTitle() + "-cellpose")
//                IJ.saveAs(refChSeg, "Tiff", outputDir.getAbsolutePath() + File.separator + imps[j].getTitle().replaceAll("/", ""))
//                refChSeg.hide()
//
//                /** Analyzing target channel (signal) == green */
//                targetCh.show()
//                IJ.run(targetCh, "Cellpose Advanced (custom model)", String.format("diameter=0 cellproba_threshold=0.0 flow_threshold=1.5 anisotropy=1.0 diam_threshold=9999999 model_path=%s model=%s nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=", greenModel.getAbsolutePath(), greenModel.getAbsolutePath()));
//                targetCh.hide()
//                targetChSeg = WindowManager.getImage(targetCh.getShortTitle() + "-cellpose")
//                IJ.saveAs(targetChSeg, "Tiff", outputDir.getAbsolutePath() + File.separator + imps[j].getTitle().replaceAll("/", ""))
//                targetChSeg.hide()

//        }

    /** Get ref signal population (signal) */
    def imgRef = ImageInt.wrap(extractCurrentStack(refChSeg));
    def populationRef = new Objects3DPopulation(imgRef);
    /** Get target signal population (signal) */
    def imgTarget = ImageInt.wrap(extractCurrentStack(targetChSeg));
    def populationTarget = new Objects3DPopulation(imgTarget);

    /** Flip target image horizontally and vertically */
    IJ.run(targetChSeg, "Flip Horizontally", "stack");
    IJ.run(targetChSeg, "Flip Vertically", "stack");

    /** Get target signal population (signal) */
    def imgTargetFlip = ImageInt.wrap(extractCurrentStack(targetChSeg));
    def populationTargetFlip = new Objects3DPopulation(imgTargetFlip);

    for (def a = 0.intValue(); a < clones.length; a++) {
        if (refChSeg.getTitle().contains(clones[a].toString())) {
            /** Get min Distances in raw image */
            def minDistPixel = new ArrayList<Double>();
            def minDistUnit = new ArrayList<Double>();
            def minDistPixelFlip = new ArrayList<Double>();
            def minDistUnitFlip = new ArrayList<Double>();

            for (def k = 0.intValue(); k < populationRef.getNbObjects(); k++) {
                def perCoordinateDistPixel = new ArrayList<Double>();
                def perCoordinateDistUnit = new ArrayList<Double>();
                def perCoordinateDistPixelFlip = new ArrayList<Double>();
                def perCoordinateDistUnitFlip = new ArrayList<Double>();
                for (def l = 0.intValue(); l < populationTarget.getNbObjects(); l++) {
                    //if (populationGreen.getObject(k).getCenterX() != populationGreen.getObject(k).getCenterX()) {
                    perCoordinateDistPixel.add(populationRef.getObject(k).distBorderPixel(populationTarget.getObject(l)))
                    perCoordinateDistUnit.add(populationRef.getObject(k).distBorderUnit(populationTarget.getObject(l)))
                    perCoordinateDistPixelFlip.add(populationRef.getObject(k).distBorderPixel(populationTargetFlip.getObject(l)))
                    perCoordinateDistUnitFlip.add(populationRef.getObject(k).distBorderUnit(populationTargetFlip.getObject(l)))
                    //}
                }
                if (perCoordinateDistPixel.size() > 1)
                    minDistPixel.add(Collections.min(perCoordinateDistPixel))
                else if (perCoordinateDistPixel.size() == 1)
                    minDistPixel.add(perCoordinateDistPixel.get(0))
                if (perCoordinateDistUnit.size() > 1)
                    minDistUnit.add(Collections.min(perCoordinateDistUnit))
                else if (perCoordinateDistUnit.size() == 1)
                    minDistUnit.add(perCoordinateDistUnit.get(0))
                if (perCoordinateDistPixelFlip.size() > 1)
                    minDistPixelFlip.add(Collections.min(perCoordinateDistPixelFlip))
                else if (perCoordinateDistPixelFlip.size() == 1)
                    minDistPixelFlip.add(perCoordinateDistPixelFlip.get(0))
                if (perCoordinateDistUnitFlip.size() > 1)
                    minDistUnitFlip.add(Collections.min(perCoordinateDistUnitFlip))
                else if (perCoordinateDistUnitFlip.size() == 1)
                    minDistUnitFlip.add(perCoordinateDistUnit.get(0))
            }
        }
    }
    /** Save arraylist as arrays to get into KolmogorovSmirnov */
    IJ.log(minDistPixel.size() + "-------------------mindistpixel")
    IJ.log(minDistUnit.size() + "-------------------mindistunit")
    IJ.log(minDistPixelFlip.size() + "-------------------mindistpixelflip")
    IJ.log(minDistUnitFlip.size() + "-------------------mindistunitflip")

    def minDistPixelAr = Stream.of(minDistPixel.stream().toArray(Double[]::new)).mapToDouble(Double::doubleValue).toArray();
    def minDistUnitAr = Stream.of(minDistUnit.stream().toArray(Double[]::new)).mapToDouble(Double::doubleValue).toArray();
    def minDistPixelFlipAr = Stream.of(minDistPixelFlip.stream().toArray(Double[]::new)).mapToDouble(Double::doubleValue).toArray();
    def minDistUnitFlipAr = Stream.of(minDistUnitFlip.stream().toArray(Double[]::new)).mapToDouble(Double::doubleValue).toArray();

    /** KolmogorovSmirnov test */
    def kst = new KolmogorovSmirnovTest()
    def pValuePixel = 0.intValue()
    def pValueUnit = 0.intValue()
    if (minDistPixel.size() != 0 && minDistPixelFlip.size() != 0)
        pValuePixel = kst.kolmogorovSmirnovStatistic(minDistPixelAr, minDistPixelFlipAr)
    else
        pValuePixel = "null"
    if (minDistUnit.size() != 0 && minDistUnitFlip.size() != 0)
        pValueUnit = kst.kolmogorovSmirnovStatistic(minDistUnitAr, minDistUnitFlipAr)
    else
        pValueUnit = "null"
    /** Iterate through table per serie for KST p-value */
    tablePerFile.incrementCounter()
    if (refIndex == 1)
        tablePerFile.setValue("Image Serie", i, refChSeg.getTitle().replaceAll("_ch_G", ""))
    else
        tablePerFile.setValue("Image Serie", i, refChSeg.getTitle().replaceAll("_ch_B", ""))
    tablePerFile.setValue("Channel Source", i, refIndex.toString())
    tablePerFile.setValue("Channel Target", i, targetIndex.toString())
    tablePerFile.setValue("Association KST p-value (Pixel)", i, pValuePixel.toString())
    tablePerFile.setValue("Association KST p-value (Unit)", i, pValueUnit.toString())


    /** Create overlap table per label  */

    //totalTable.setPrecision( 10 );

    def refLabelIDs = LabelImages.findAllLabels(refChSeg)

    for (def m = 0.intValue(); m < refLabelIDs.length; m++) {

        tableOverlap.incrementCounter()
        // Image Serie
        if (refIndex == 1)
            tableOverlap.setValue("Image Serie", m, refChSeg.getTitle().replaceAll("_ch_G", ""))
        else
            tableOverlap.setValue("Image Serie", m, refChSeg.getTitle().replaceAll("_ch_B", ""))
        // Channel Source
        tableOverlap.setValue("Channel Source", m, refIndex.toString())
        // Channel Target
        tableOverlap.setValue("Channel Target", m, targetIndex.toString())
        // Label ID
        tableOverlap.setValue("Label ID (Source)", m, refLabelIDs[m].toString())
        // Overlap
        tableOverlap.setValue("Target Overlap", m, LabelImages.getTargetOverlapPerLabel(refChSeg, targetChSeg).getStringValue("TargetOverlap", m));
        // Jaccard index
        tableOverlap.setValue("Jaccard Index", m, LabelImages.getJaccardIndexPerLabel(refChSeg, targetChSeg).getStringValue("JaccardIndex", m));
        // Dice coefficient
        tableOverlap.setValue("Dice Coefficient", m, LabelImages.getDiceCoefficientPerLabel(refChSeg, targetChSeg).getStringValue("DiceCoefficient", m));
        // Volume similarity
        tableOverlap.setValue("Volume Similarity", m, LabelImages.getVolumeSimilarityPerLabel(refChSeg, targetChSeg).getStringValue("VolumeSimilarity", m));
        if (LabelImages.getFalseNegativeErrorPerLabel(refChSeg, targetChSeg) != null)
        // False negative error
            tableOverlap.setValue("False Negative Error", m, LabelImages.getFalseNegativeErrorPerLabel(refChSeg, targetChSeg).getStringValue("FalseNegativeError", m));
        else
            tableOverlap.setValue("False Negative Error", m, "null");
        if (LabelImages.getFalsePositiveErrorPerLabel(refChSeg, targetChSeg) != null)
        // False positive error
            tableOverlap.setValue("False Positive Error", m, LabelImages.getFalsePositiveErrorPerLabel(refChSeg, targetChSeg).getStringValue("FalsePositiveError", m));
        else
            tableOverlap.setValue("False Positive Error", m, "null");
    }
    /** Create overlap table per serie  */
    totalTable.incrementCounter();
    // Image Serie
    if (refIndex == 1)
        totalTable.setValue("Image Serie", i, refChSeg.getTitle().replaceAll("_ch_G", ""));
    else
        totalTable.setValue("Image Serie", i, refChSeg.getTitle().replaceAll("_ch_B", ""));
    // Channel Source
    totalTable.setValue("Channel Source", i, refIndex.toString());

    // Channel Target
    totalTable.setValue("Channel Target", i, targetIndex.toString());

    // Overlap
    totalTable.setValue("Total Overlap", i, LabelImages.getTotalOverlap(refChSeg, targetChSeg));

    // Jaccard index
    totalTable.setValue("Jaccard Index", i, LabelImages.getJaccardIndex(refChSeg, targetChSeg));

    // Dice coefficient
    totalTable.setValue("Dice Coefficient", i, LabelImages.getDiceCoefficient(refChSeg, targetChSeg));

    // Volume similarity
    totalTable.setValue("Volume Similarity", i, LabelImages.getVolumeSimilarity(refChSeg, targetChSeg));

    // False negative error
    totalTable.setValue("False Negative Error", i, LabelImages.getFalseNegativeError(refChSeg, targetChSeg));

    // False positive error
    totalTable.setValue("False Positive Error", i, LabelImages.getFalsePositiveError(refChSeg, targetChSeg));
//        /** Measure morphology on 3D labels, get ResultsTable */
//        def morphoTableRef = new inra.ijpb.plugins.AnalyzeRegions3D().process(refChSeg)
//        def morphoTableTarget = new inra.ijpb.plugins.AnalyzeRegions3D().process(targetChSeg)
//        /** Save morphology table  per serie  */
//        morphoTableRef.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[i].getName() + "_morphology_measurements" + ".csv")
//        morphoTableTarget.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesTarget[i].getName() + "_morphology_measurements" + ".csv")

    //}

    /** Save overlap table per label per serie per file */
    // set 6 decimal places in the displayed results
    tableOverlap.setPrecision(6);
    if (refIndex == 1)
        tableOverlap.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[i].getName().replaceAll("_ch_G", "") + "_individual_labels_overlap_measurements" + ".csv")
    else
        tableOverlap.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[i].getName().replaceAll("_ch_B", "") + "_individual_labels_overlap_measurements" + ".csv")

    //}


}
/** Save KST table per serie per file */
// set 6 decimal places in the displayed results
tablePerFile.setPrecision(6);
if (refIndex == 1)
    tablePerFile.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[0].getName() + "_all_labels_KST_pvalue" + ".csv")
else
    tablePerFile.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[0].getName() + "_all_labels_KST_pvalue" + ".csv")

//set 6 decimal places in the displayed results
totalTable.setPrecision(6);
/** Save overlap table  per serie per file */
if (refIndex == 1)
    totalTable.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[0].getName() + "_all_labels_overlap_measurements" + ".csv")
else
    totalTable.saveAs(outputDir.getAbsolutePath() + File.separator + listOfFilesRef[0].getName() + "_all_labels_overlap_measurements" + ".csv")

IJ.log("Done!!!")

ImagePlus extractCurrentStack(ImagePlus plus) {
    // check dimensions
    int[] dims = plus.getDimensions();//XYCZT
    int channel = plus.getChannel();
    int frame = plus.getFrame();
    ImagePlus stack;
    // crop actual frame
    if ((dims[2] > 1) || (dims[4] > 1)) {
        IJ.log("hyperstack found, extracting current channel " + channel + " and frame " + frame);
        def duplicator = new Duplicator();
        stack = duplicator.run(plus, channel, channel, 1, dims[3], frame, frame);
    } else stack = plus.duplicate();

    return stack;
}

