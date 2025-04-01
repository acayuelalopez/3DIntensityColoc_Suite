import ij.IJ
import ij.ImagePlus
import ij.gui.Roi
import ij.gui.ShapeRoi
import ij.measure.ResultsTable
import ij.plugin.ChannelSplitter
import ij.plugin.Duplicator
import ij.plugin.RGBStackMerge
import org.apache.commons.math3.stat.inference.TTest
import mcib3d.geom.Object3D
import mcib3d.geom.Objects3DPopulation
import mcib3d.image3d.ImageInt

// INPUT UI
//
//#@File(label = "Input Segmentation Files Directory", style = "directory") inputFilesSeg
//#@File(label = " Input Raw Files Directory", style = "directory") inputFilesRaw
//#@File(label = "output directory", style = "directory") outputDir
//#@Integer(label = "Nuclei channel", value = 2) nucleiChannel
//#@Integer(label = "Telomere channel", value = 0) telomereChannel
//#@Integer(label = "Marker channel", value = 1) markerChannel


// IDE
//
//
def inputFilesTrf1 = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/telomere")
def inputFilesNuclei = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/nuclei")
def inputFilesRaw = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output/images")
def outputDir = new File("/mnt/imgserver/CONFOCAL/IA/Projects/2025/2025_2_7_jfiel/output")
def nucleiChannel = 0.intValue()
def trf1Channel = 1.intValue()

//def headless = true;
//new ImageJ().setVisible(true);

IJ.log("-Parameters selected: ")
IJ.log("    -Input Seg Files Dir: " + inputFilesTrf1)
IJ.log("    -Input Raw Files Dir: " + inputFilesRaw)
IJ.log("    -output Dir: " + outputDir)
IJ.log("    -Nuclei Channel: " + nucleiChannel)
IJ.log("    -Telomere Channel: " + trf1Channel)
IJ.log("                                                           ");
/** Get files (images) from input directory */
def listOfFiles = inputFilesRaw.listFiles();
def wtRt = new ResultsTable()
def koRt = new ResultsTable()
def pValueRt = new ResultsTable()


Map<String, Integer> patternCount = new HashMap<>();

for (File file : listOfFiles) {
    String pattern = extractPattern(file.getName());
    patternCount.put(pattern, patternCount.getOrDefault(pattern, 0) + 1);
}

List<String> repeatedPatterns = new ArrayList<>();
for (Map.Entry<String, Integer> entry : patternCount.entrySet()) {
    if (entry.getValue() > 1) {
        repeatedPatterns.add(entry.getKey());
    }
}


def counterWt = 0.intValue()
def counterKo = 0.intValue()
for (def x = 0; x < repeatedPatterns.size(); x++) {

    //WT
/////Trf1
    def wtTrf1N_All = new ArrayList<Double>()
    def wtTrf1MeanInt_All = new ArrayList<Double>()
    def wtTrf1SumInt_All = new ArrayList<Double>()
    def wtTrf1StdInt_All = new ArrayList<Double>()

//KO
//////Trf1
    def koTrf1N_All = new ArrayList<Double>()
    def koTrf1MeanInt_All = new ArrayList<Double>()
    def koTrf1SumInt_All = new ArrayList<Double>()
    def koTrf1StdInt_All = new ArrayList<Double>()



    for (def i = 0; i < listOfFiles.length; i++) {
        //WT
        def wtNucleiN = new ArrayList<Double>()
//Area
        def wtNucleusArea = new ArrayList<Double>()
//TRF1
        def wtTrf1N = new ArrayList<Double>()
        def wtTrf1MeanInt = new ArrayList<Double>()
        def wtTrf1SumInt = new ArrayList<Double>()
        def wtTrf1StdInt = new ArrayList<Double>()


        //KO
        def koNucleiN = new ArrayList<Double>()
//Area
        def koNucleusArea = new ArrayList<Double>()
//TRF1
        def koTrf1N = new ArrayList<Double>()
        def koTrf1MeanInt = new ArrayList<Double>()
        def koTrf1SumInt = new ArrayList<Double>()
        def koTrf1StdInt = new ArrayList<Double>()



        def counter = 0.intValue()
        def tablePerImage = new ResultsTable();
        /** Create image for each file in the input directory */
        def imps = new ImagePlus(inputFilesRaw.getAbsolutePath() + File.separator + listOfFiles[i].getName())
        def cal = imps.getCalibration()
        IJ.log(imps.getTitle())
        /** Split channels */
        def channels = ChannelSplitter.split(imps)

        /** Get trf1 channel */
        def labelTrf1 = new ImagePlus(inputFilesTrf1.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
        def chTrf1ToMeasure = channels[trf1Channel.intValue()]

        /** Get nuclei channel */
        def labelNuclei = new ImagePlus(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))
        labelNuclei.setCalibration(cal)
        def chNucleiToMeasure = channels[nucleiChannel.intValue()]
        IJ.log(inputFilesNuclei.getAbsolutePath() + File.separator + listOfFiles[i].getName().replaceAll(".tif", "_cp_masks.tif"))

        // Get Nuclei objects population
        def imgNuclei = ImageInt.wrap(extractCurrentStack(labelNuclei));
        def populationNuclei = new Objects3DPopulation(imgNuclei);
        // Get Nuclei signal
        def signalNuclei = ImageInt.wrap(extractCurrentStack(chNucleiToMeasure));

        // Get Trf1 objects population
        def imgTrf1 = ImageInt.wrap(extractCurrentStack(labelTrf1));
        def populationTrf1 = new Objects3DPopulation(imgTrf1);
        // Get Trf1 signal
        def signalTrf1 = ImageInt.wrap(extractCurrentStack(chTrf1ToMeasure));

        IJ.saveAs(RGBStackMerge.mergeChannels(new ImagePlus[]{labelNuclei, chNucleiToMeasure, labelTrf1, chTrf1ToMeasure}, false), "Tiff", outputDir.getAbsolutePath() + File.separator + "merge" + File.separator + listOfFiles[i].getName())

        if (listOfFiles[i].getName().contains("WT")) {


            if (listOfFiles[i].getName().contains(repeatedPatterns.get(x))) {
                wtNucleiN.add(populationNuclei.getNbObjects().doubleValue())
                IJ.log("pasaaaaaaaaaaaa 0000000")
                for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
                    def zRange = (populationNuclei.getObject(j).zmax - populationNuclei.getObject(j).zmin)
                    if (zRange > 2) {
                        def counterTrf1 = 0.intValue()
                        def meanIntTrf1 = new ArrayList<Double>()
                        def sumIntTrf1 = new ArrayList<Double>()
                        def stdIntTrf1 = new ArrayList<Double>()


                        wtNucleusArea.add(populationNuclei.getObject(j).volumeUnit)
                        for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                            if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint())) {
                                counterTrf1++
                                meanIntTrf1.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                                sumIntTrf1.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                                stdIntTrf1.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                                //wtTrf1MeanInt_All.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                                //wtTrf1SumInt_All.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                                //wtTrf1StdInt_All.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                            }
                        }


                        wtTrf1N.add(counterTrf1.doubleValue())
                        wtTrf1MeanInt.add(meanIntTrf1.stream()
                                .mapToDouble(d -> d)
                                .average()
                                .orElse(0.0))

                        wtTrf1SumInt.add(sumIntTrf1.stream()
                                .mapToDouble(d -> d)
                                .sum())

                        wtTrf1StdInt.add(stdIntTrf1.stream()
                                .mapToDouble(d -> d)
                                .average()
                                .orElse(0.0))


                    }
                }

                wtTrf1N_All.add(wtTrf1N.stream()
                        .mapToDouble(d -> d)
                        .sum())

                wtTrf1MeanInt_All.add(wtTrf1MeanInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))
                wtTrf1SumInt_All.add(wtTrf1SumInt.stream()
                        .mapToDouble(d -> d)
                        .sum())
                wtTrf1StdInt_All.add(wtTrf1StdInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))



            }

        }

        if (listOfFiles[i].getName().contains("KO")) {


            if (listOfFiles[i].getName().contains(repeatedPatterns.get(x))) {

                koNucleiN.add(populationNuclei.getNbObjects().doubleValue())
                IJ.log("pasaaaaaaaaaaaa 0000000")
                for (def j = 0.intValue(); j < populationNuclei.getNbObjects(); j++) {
                    def zRange = (populationNuclei.getObject(j).zmax - populationNuclei.getObject(j).zmin)
                    if (zRange > 2) {
                        def counterTrf1 = 0.intValue()
                        def meanIntTrf1 = new ArrayList<Double>()
                        def sumIntTrf1 = new ArrayList<Double>()
                        def stdIntTrf1 = new ArrayList<Double>()


                        koNucleusArea.add(populationNuclei.getObject(j).volumeUnit)
                        for (def k = 0.intValue(); k < populationTrf1.getNbObjects(); k++) {
                            if (populationNuclei.getObject(j).inside(populationTrf1.getObject(k).getCenterAsPoint())) {
                                counterTrf1++
                                meanIntTrf1.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                                sumIntTrf1.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                                stdIntTrf1.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                                //wtTrf1MeanInt_All.add(populationTrf1.getObject(k).getPixMeanValue(signalTrf1))
                                //wtTrf1SumInt_All.add(populationTrf1.getObject(k).getIntegratedDensity(signalTrf1))
                                //wtTrf1StdInt_All.add(populationTrf1.getObject(k).getPixStdDevValue(signalTrf1))
                            }
                        }


                        koTrf1N.add(counterTrf1.doubleValue())
                        koTrf1MeanInt.add(meanIntTrf1.stream()
                                .mapToDouble(d -> d)
                                .average()
                                .orElse(0.0))

                        koTrf1SumInt.add(sumIntTrf1.stream()
                                .mapToDouble(d -> d)
                                .sum())

                        koTrf1StdInt.add(stdIntTrf1.stream()
                                .mapToDouble(d -> d)
                                .average()
                                .orElse(0.0))


                    }
                }

                koTrf1N_All.add(koTrf1N.stream()
                        .mapToDouble(d -> d).
                        sum())

                koTrf1MeanInt_All.add(koTrf1MeanInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))
                koTrf1SumInt_All.add(koTrf1SumInt.stream()
                        .mapToDouble(d -> d).
                        sum())
                koTrf1StdInt_All.add(koTrf1StdInt.stream()
                        .mapToDouble(d -> d)
                        .average()
                        .orElse(0.0))


            }

        }


    }
    wtRt.incrementCounter()
    counterWt++
    wtRt.setValue("Animal", counterWt, repeatedPatterns.get(x))
    wtRt.setValue("N of Telomeres per Nucleus", counterWt, wtTrf1N_All.stream()
            .mapToDouble(d -> d)
            .sum());
    ;
    wtRt.setValue("Mean Intensity of Telomeres per Nucleus", counterWt, wtTrf1MeanInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0));

    wtRt.setValue("Sum Intensity of Telomeres per Nucleus", counterWt, wtTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .sum());


    wtRt.setValue("Mean of Sum Intensity of Telomeres per Nucleus", counterWt, wtTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0))
    wtRt.setValue("Std of Sum Intensity of Telomeres per Nucleus", counterWt, std(wtTrf1SumInt_All, wtTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0)))
    wtRt.setValue("SEM of Sum Intensity of Telomeres per Nucleus", counterWt, std(wtTrf1SumInt_All, wtTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0)) / Math.sqrt(wtTrf1SumInt_All.size()))
    koRt.incrementCounter()
    counterKo++
    koRt.setValue("Animal", counterKo, repeatedPatterns.get(x))
    koRt.setValue("N of Telomeres per Nucleus", counterKo, koTrf1N_All.stream()
            .mapToDouble(d -> d)
            .sum());

    koRt.setValue("Mean Intensity of Telomeres per Nucleus", counterKo, koTrf1MeanInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0));
    koRt.setValue("Sum Intensity of Telomeres per Nucleus", counterKo, koTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .sum());


    koRt.setValue("Mean of Sum Intensity of Telomeres per Nucleus", counterKo, koTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0));
    koRt.setValue("Std of Sum Intensity of Telomeres per Nucleus", counterKo, std(koTrf1SumInt_All, koTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0)))
    koRt.setValue("SEM of Sum Intensity of Telomeres per Nucleus", counterKo, std(koTrf1SumInt_All, koTrf1SumInt_All.stream()
            .mapToDouble(d -> d)
            .average()
            .orElse(0.0)) / Math.sqrt(koTrf1SumInt_All.size()))
}

/*def tTest = new TTest()
IJ.log(wtTrf1N_All.size() + "------" + wtTrf1N_All + "---" + koTrf1N_All.size() + "----" + koTrf1N_All)
pValueRt.setValue("p-value N of TRF1 t-Test", 0, tTest.tTest(wtTrf1N_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1N_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value N of Colocs t-Test", 0, tTest.tTest(wtColocsN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koColocsN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Mean Intensity TRF1 t-Test", 0, new BigDecimal(tTest.tTest(wtTrf1MeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1MeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray())).toString())
pValueRt.setValue("p-value Sum Intensity TRF1 t-Test", 0, tTest.tTest(wtTrf1SumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1SumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Std Intensity TRF1 t-Test", 0, tTest.tTest(wtTrf1StdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koTrf1StdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value N of H2AX t-Test", 0, tTest.tTest(wtH2axN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axN_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Mean Intensity H2AX t-Test", 0, tTest.tTest(wtH2axMeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axMeanInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Sum Intensity H2AX t-Test", 0, tTest.tTest(wtH2axSumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axSumInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.setValue("p-value Std Intensity H2AX t-Test", 0, tTest.tTest(wtH2axStdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray(), koH2axStdInt_All.stream()
        .mapToDouble(Double::doubleValue)
        .toArray()).toString())
pValueRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3D_pvalue_tTest_WT_KO_per_Clon.csv")*/
wtRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_WT_per_Clon.csv")
koRt.saveAs(outputDir.getAbsolutePath() + File.separator + "csv" + File.separator + "3DAnalysis_KO_per_Clon.csv")

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

static double std(ArrayList<Double> table, double mean) {
    // Step 1:
    double meanDef = mean
    double temp = 0;

    for (int i = 0; i < table.size(); i++) {
        int val = table.get(i);

        // Step 2:
        double squrDiffToMean = Math.pow(val - meanDef, 2);

        // Step 3:
        temp += squrDiffToMean;
    }

    // Step 4:
    double meanOfDiffs = (double) temp / (double) (table.size());

    // Step 5:
    return Math.sqrt(meanOfDiffs);
}

private static String extractPattern(String fileName) {
    int startIndex = fileName.indexOf("KO");
    if (startIndex == -1) {
        startIndex = fileName.indexOf("WT");
    }
    int endIndex = fileName.indexOf("Series", startIndex) + "Series".length();
    return fileName.substring(startIndex, endIndex);
}
