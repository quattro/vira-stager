package edu.gsu.cs.vira;

import ch.ethz.bsse.indelfixer.minimal.Start;
import net.sf.samtools.SAMFileReader;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import qure.AmpliconSet;
import qure.ReadSet;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class Main {

    @Option(name = "-i")
    private String input;

    @Option(name = "-g")
    private String reference;

    @Option(name = "-k")
    private int numHap = 20;

    @Option(name = "-o", usage = "Path to the output directory (default: current directory)")
    private String output = "./";

    @Option(name = "-n", usage = "Number of iterations for interval finding")
    private int iterations = 1000;

    @Option(name = "-p", usage = "Number of processors to use")
    private int threads = Runtime.getRuntime().availableProcessors();

    private static final String IDF_OUTPUT = "reads.sam";

    public static void main(String[] args) {
        new Main().doMain(args);
    }

    private void doMain(String[] args) {
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);

            if (this.input == null && this.reference == null) {
                throw new CmdLineException("");
            }
            if (this.output == null) {
                this.output = "";
            }
            if (this.iterations <= 0 || this.threads <= 0) {
                throw new CmdLineException("");
            }

            Start idf = new Start();
            String[] idfArgs;
            if (this.output.isEmpty()) {
                idfArgs = new String[]{"-i", this.input, "-g", this.reference,  "-rmDel", "-454"};
            } else {
                idfArgs = new String[]{"-i", this.input, "-g", this.reference, "-o", this.output,  "-rmDel", "-454"};
            }
            idf.doMain(idfArgs);
            String idfOutputPath = this.output.isEmpty() ? IDF_OUTPUT : this.output + File.separator + IDF_OUTPUT;

            File samFile = new File(idfOutputPath);
            File refFile = new File(this.reference);
            DNASequence reference = FastaReaderHelper.readFastaDNASequence(refFile, false).values().iterator().next();
            SAMFileReader samReader = new SAMFileReader(samFile);

            // Use the QuRe random-algorithm to find good intervals on IDF alignment
            ReadSet readSet = new ReadSet();
            readSet.setPopulation(samReader, reference);
            readSet.buildDictionary(9);
            readSet.estimateAmpliconsParallel(this.iterations, this.threads);
            AmpliconSet intervals = readSet.ampliconSet;

            samReader.close();

            double[] starts = intervals.getStarts();
            double[] stops = intervals.getStops();

            // Use kGEM on each interval to correct reads.
            for (int i = 0; i < starts.length; i++) {
                int start = (int) starts[i];
                int stop = (int) stops[i];

                String ampDirPathName = this.output.isEmpty() ? "amplicon" + i : this.output + File.separator + "amplicon" + i;
                Path ampDirPath = Paths.get(ampDirPathName);
                if (!Files.exists(ampDirPath))
                    Files.createDirectory(ampDirPath);

                Path ampFASTAFilePath = Paths.get(ampDirPath.toString() + File.separator + "reads.fa");
                if (Files.exists(ampFASTAFilePath))
                    Files.delete(ampFASTAFilePath);

                String[] erifArgs = new String[]{
                        "-sam", idfOutputPath,
                        "-g", this.reference,
                        "-r", String.valueOf(start) + "-" + String.valueOf(stop),
                        "-o", ampDirPath.toString() + File.separator};
                edu.gsu.cs.align.exec.Main.main(erifArgs);

                String correctedFASTAPathName = ampDirPathName + File.separator + "aligned_reads.fas";
                String[] kgemArgs = new String[]{correctedFASTAPathName, Integer.toString(numHap),
                        "-o", ampDirPathName + File.separator + "corrected.fa", "-r"};
                edu.gsu.cs.kgem.exec.Main.main(kgemArgs);
            }


        } catch (CmdLineException e) {
           parser.printUsage(System.err);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);

        }
        System.exit(0);
    }

    static void printHelp() {
        System.out.println("");
    }
}