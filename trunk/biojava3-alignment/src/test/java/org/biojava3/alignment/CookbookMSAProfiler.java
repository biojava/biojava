package org.biojava3.alignment;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava3.alignment.Alignments.ProfileProfileAlignerType;
import org.biojava3.alignment.template.GapPenalty;
import org.biojava3.alignment.template.PairwiseSequenceScorer;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.io.FastaReaderHelper;
import org.biojava3.core.util.ConcurrencyTools;

public class CookbookMSAProfiler {

    private static class Profiler {

        private long maxMemoryUsed, timeCheckpoint;
        private final long timeStart;

        private Profiler() {
            maxMemoryUsed = Runtime.getRuntime().totalMemory();
            timeStart = timeCheckpoint = System.nanoTime();
        }

        private long getMaxMemoryUsed() {
            return maxMemoryUsed = Math.max(maxMemoryUsed, Runtime.getRuntime().totalMemory());
        }

        private long getTimeSinceCheckpoint() {
            return System.nanoTime() - timeCheckpoint;
        }

        private long getTimeSinceStart() {
            return System.nanoTime() - timeStart;
        }

        private void setCheckpoint() {
            maxMemoryUsed = Math.max(maxMemoryUsed, Runtime.getRuntime().totalMemory());
            timeCheckpoint = System.nanoTime();
        }

    }

    public static void main(String[] args) throws Exception {

        if (args.length < 1) {
            System.err.println("The first argument must be a fasta file of protein sequences.");
            return;
        }

        // ConcurrencyTools.setThreadPoolSingle();

        PrintStream fout = new PrintStream("msa.txt");
        Profiler profiler = new Profiler();

        System.out.printf("Loading sequences from %s... ", args[0]);
        List<ProteinSequence> list = new ArrayList<ProteinSequence>();
        list.addAll(FastaReaderHelper.readFastaProteinSequence(new File(args[0])).values());
        if (args.length > 1 && Integer.parseInt(args[1]) < list.size()) {
            System.out.printf("%s/%d", args[1], list.size());
            list = list.subList(0, Integer.parseInt(args[1]));
        } else {
            System.out.printf("%d", list.size());
        }
        System.out.printf(" sequences in %d ms using %d kB%n%n", profiler.getTimeSinceCheckpoint()/1000000,
                profiler.getMaxMemoryUsed()/1024);

        profiler.setCheckpoint();

        System.out.print("Stage 1: pairwise similarity calculation... ");
        GapPenalty gaps = new SimpleGapPenalty();
        SubstitutionMatrix<AminoAcidCompound> blosum62 = new SimpleSubstitutionMatrix<AminoAcidCompound>();
        List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> scorers = Alignments.getAllPairsScorers(list,
                PairwiseSequenceScorerType.GLOBAL_IDENTITIES, gaps, blosum62);
        Alignments.runPairwiseScorers(scorers);
        System.out.printf("%d scores in %d ms using %d kB%n%n", scorers.size(),
                profiler.getTimeSinceCheckpoint()/1000000, profiler.getMaxMemoryUsed()/1024);

        profiler.setCheckpoint();

        System.out.print("Stage 2: hierarchical clustering into a guide tree... ");
        GuideTree<ProteinSequence, AminoAcidCompound> tree = new GuideTree<ProteinSequence, AminoAcidCompound>(list,
                scorers);
        scorers = null;
        System.out.printf("%d ms using %d kB%n%n%s%n%n", profiler.getTimeSinceCheckpoint()/1000000,
                profiler.getMaxMemoryUsed()/1024, tree);

        profiler.setCheckpoint();

        System.out.print("Stage 3: progressive alignment... ");
        Profile<ProteinSequence, AminoAcidCompound> msa = Alignments.getProgressiveAlignment(tree,
                ProfileProfileAlignerType.GLOBAL, gaps, blosum62);
        System.out.printf("%d profile-profile alignments in %d ms using %d kB%n%n", list.size() - 1,
                profiler.getTimeSinceCheckpoint()/1000000, profiler.getMaxMemoryUsed()/1024);
        fout.print(msa);
        fout.close();

        ConcurrencyTools.shutdown();

        System.out.printf("Total time: %d ms%nMemory use: %d kB%n", profiler.getTimeSinceStart()/1000000,
                profiler.getMaxMemoryUsed()/1024);

    }

}
