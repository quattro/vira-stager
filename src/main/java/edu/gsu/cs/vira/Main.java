package edu.gsu.cs.vira;

import ch.ethz.bsse.indelfixer.minimal.Start;
import net.sf.picard.sam.BuildBamIndex;
import net.sf.picard.sam.SamFormatConverter;
import net.sf.picard.sam.SortSam;
import net.sf.samtools.SAMFileReader;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.FastaSequenceParser;
import org.biojava3.core.sequence.io.FileProxyDNASequenceCreator;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import qure.AmpliconSet;
import qure.ReadSet;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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


    public static void main(String[] args) {
        new Main().doMain(args);
    }

    private void doMain(String[] args) {
        CmdLineParser parser = new CmdLineParser(this);
        BufferedWriter intervalWriter  = null;
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


            File samFile = new File(this.input);
            File refFile = new File(this.reference);
            // Use this instead of the helper. This way we can parse ambiguous base-codes
            FastaReader<DNASequence, NucleotideCompound > fastaProxyReader =
                    new FastaReader<>(
                            refFile,
                            new GenericFastaHeaderParser<DNASequence, NucleotideCompound>(),
                            new FileProxyDNASequenceCreator(
                                    refFile,
                                    AmbiguityDNACompoundSet.getDNACompoundSet(),
                                    new FastaSequenceParser()
                            )
                    );
            DNASequence reference = null;
            for (DNASequence seq : fastaProxyReader.process().values()) {
                reference = seq;
                break;
            }

            String bamName = this.input.replace("sam", "bam");
            String sortedName = bamName.replace("bam", "sorted.bam");

            new SamFormatConverter().instanceMain(new String[]{"Input=" + this.input,
                    "Output=" + bamName, "VERBOSITY=ERROR", "QUIET=TRUE"});
            new SortSam().instanceMain(new String[]{"Input=" + bamName, "Output=" + sortedName,
                    "SORT_ORDER=coordinate", "VERBOSITY=ERROR", "QUIET=TRUE"});
            new BuildBamIndex().instanceMain(new String[]{"Input="+sortedName,
                    "Output="+sortedName+".bai", "VERBOSITY=ERROR", "QUIET=TRUE"});

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

            // write out first. makes it easier to debug.
            intervalWriter = new BufferedWriter(new FileWriter("intervals.txt"));
            for (int i = 0; i < starts.length; i++) {
                int start = (int) starts[i];
                int stop = (int) stops[i];

                intervalWriter.write(start + "," + stop);
                intervalWriter.newLine();
            }
            intervalWriter.close();

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

                String inter = "CONSENSUS:" + String.valueOf(start) + "-" + String.valueOf(stop);
                ProcessBuilder pb = new ProcessBuilder("b2w", sortedName, this.reference, inter,
                        "-w", String.valueOf(stop - start), "-i", "0", "-m", "20", "-x", "10000");
                Process proc = pb.start();
                proc.waitFor();

                //move the output
                String correctedFASTAPathName = ampDirPathName + File.separator + "aligned_reads.fas";
                new File("w-" + inter.replace(':', '-') + ".reads.fas").renameTo(
                                    new File(correctedFASTAPathName));

                String[] kgemArgs = new String[]{correctedFASTAPathName, Integer.toString(numHap),
                        "-o", ampDirPathName + File.separator + "corrected.fa", "-r",
                        "-t", "0", "-d", "1"};
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
}