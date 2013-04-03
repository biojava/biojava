/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.phylo;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;
import org.biojava3.core.sequence.MultipleSequenceAlignment;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;


import org.forester.io.writers.PhylogenyWriter;
import org.forester.phylogeny.Phylogeny;
import org.forester.evoinference.matrix.distance.BasicSymmetricalDistanceMatrix;
import org.forester.evoinference.matrix.distance.DistanceMatrix;
import org.forester.evoinference.distance.NeighborJoining;

/**
 * Tree constructor uses the forrester tree library to build phylogenetic trees using neighbor joining algorithm. The distance matrix
 * is calculated using code from JalView.
 *
 * @author Scooter Willis
 */
public class TreeConstructor<C extends AbstractSequence<D>, D extends Compound> extends Thread {

    TreeType treeType;
    TreeConstructionAlgorithm treeConstructionAlgorithm;
    NJTreeProgressListener treeProgessListener;
    MultipleSequenceAlignment<C, D> multipleSequenceAlignment = new MultipleSequenceAlignment<C, D>();

    public TreeConstructor(MultipleSequenceAlignment<C, D> multipleSequenceAlignment, TreeType _treeType, TreeConstructionAlgorithm _treeConstructionAlgorithm, NJTreeProgressListener _treeProgessListener) {
        treeType = _treeType;
        treeConstructionAlgorithm = _treeConstructionAlgorithm;
        treeProgessListener = _treeProgessListener;
        this.multipleSequenceAlignment = multipleSequenceAlignment;

    }

    public TreeConstructor(BasicSymmetricalDistanceMatrix _matrix, TreeType _treeType, TreeConstructionAlgorithm _treeConstructionAlgorithm, NJTreeProgressListener _treeProgessListener) {
        matrix = _matrix;
        copyDistanceMatrix = CheckTreeAccuracy.copyMatrix(matrix);
        treeType = _treeType;
        treeConstructionAlgorithm = _treeConstructionAlgorithm;
        treeProgessListener = _treeProgessListener;


    }

    public void outputPhylipDistances(String fileName) throws Exception {
        DistanceMatrix distances = getDistanceMatrix();
        if (distances == null) {
            throw new Exception("distance matrix has not been calculated. Requires process() method to be called first");
        }
        FileOutputStream fo = new FileOutputStream(fileName);
        PrintStream dos = new PrintStream(fo);
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(5);
        df.setMinimumFractionDigits(5);
        for (int row = 0; row < distances.getSize(); row++) {
            dos.print(distances.getIdentifier(row));
            for (int col = 0; col < distances.getSize(); col++) {
                dos.print(" " + df.format(distances.getValue(col, row)));
            }
            dos.println();
        }
        dos.close();
        fo.close();
    }

    private double[][] calculateDistanceMatrix(MultipleSequenceAlignment<C, D> multipleSequenceAlignment, TreeConstructionAlgorithm tca) {
        updateProgress("Determing Distances", 0);
        int numberOfSequences = multipleSequenceAlignment.getSize();
        String[] sequenceString = new String[numberOfSequences];
        for (int i = 0; i < multipleSequenceAlignment.getSize(); i++) {
            sequenceString[i] = multipleSequenceAlignment.getAlignedSequence(i).getSequenceAsString();

        }


        double[][] distance = new double[numberOfSequences][numberOfSequences];

        int totalloopcount = (numberOfSequences / 2) * (numberOfSequences + 1);

        if (tca == TreeConstructionAlgorithm.PID) {
            int loopcount = 0;
            for (int i = 0; i < (numberOfSequences - 1); i++) {
                updateProgress("Determining Distances", (loopcount * 100) / totalloopcount);
                for (int j = i; j < numberOfSequences; j++) {
                    loopcount++;
                    if (j == i) {
                        distance[i][i] = 0;
                    } else {
                        distance[i][j] = 100 - Comparison.PID(sequenceString[i], sequenceString[j]);

                        distance[j][i] = distance[i][j];
                    }
                }
            }
        } else {
            // Pairwise substitution score (with no gap penalties)
            ScoreMatrix pwmatrix = ResidueProperties.getScoreMatrix(treeConstructionAlgorithm.name());
            if (pwmatrix == null) {
                pwmatrix = ResidueProperties.getScoreMatrix(treeConstructionAlgorithm.BLOSUM62.name());
            }
            int maxscore = 0;
            int end = sequenceString[0].length();
            int loopcount = 0;
            for (int i = 0; i < (numberOfSequences - 1); i++) {
                updateProgress("Determining Distances", (loopcount * 100) / totalloopcount);
                for (int j = i; j < numberOfSequences; j++) {
                    int score = 0;
                    loopcount++;
                    for (int k = 0; k < end; k++) {
                        try {
                            score += pwmatrix.getPairwiseScore(sequenceString[i].charAt(k), sequenceString[j].charAt(k));
                        } catch (Exception ex) {
                            System.err.println("err creating BLOSUM62 tree");
                            ex.printStackTrace();
                        }
                    }

                    distance[i][j] = (float) score;

                    if (score > maxscore) {
                        maxscore = score;
                    }
                }
            }

            for (int i = 0; i < (numberOfSequences - 1); i++) {
                for (int j = i; j < numberOfSequences; j++) {
                    distance[i][j] = (float) maxscore - distance[i][j];
                    distance[j][i] = distance[i][j];
                }
            }

        }
        updateProgress("Determining Distances", 100);

        return distance;
    }

    public DistanceMatrix getDistanceMatrix() {
        return copyDistanceMatrix;
    }

    public void cancel() {
        //    if (njtree != null) {
        //        njtree.cancel();
        //    }
    }
    boolean verbose = false;
    Phylogeny p = null;
    BasicSymmetricalDistanceMatrix matrix = null;
    DistanceMatrix copyDistanceMatrix = null;

    public void process() throws Exception {


        if (matrix == null) {
            double[][] distances = calculateDistanceMatrix(multipleSequenceAlignment, treeConstructionAlgorithm);
            matrix = new BasicSymmetricalDistanceMatrix(multipleSequenceAlignment.getSize());
            for (int i = 0; i < matrix.getSize(); i++) {
                matrix.setIdentifier(i, multipleSequenceAlignment.getAlignedSequence(i).getAccession().getID());
            }
            for (int col = 0; col < matrix.getSize(); col++) {
                for (int row = 0; row < matrix.getSize(); row++) {
                    matrix.setValue(col, row, distances[col][row]);

                }
            }
            copyDistanceMatrix = CheckTreeAccuracy.copyMatrix(matrix);
        }

        final List<Phylogeny> ps = new ArrayList<Phylogeny>();
        final NeighborJoining nj = NeighborJoining.createInstance(verbose);

        ps.add(nj.execute(matrix));
        p = ps.get(0);

    }

    //   public void getTreeAccuracy(){
    //   CheckTreeAccuracy checkTreeAccuracy = new CheckTreeAccuracy();
    //   checkTreeAccuracy.process(p,distanceMatrix );
    //   }
    @Override
    public void run() {
        try {
            process();
        } catch (Exception e) {
            e.printStackTrace();

        }
    }

    public String getNewickString(boolean simpleNewick, boolean writeDistanceToParent) throws Exception {
        final PhylogenyWriter w = new PhylogenyWriter();
        StringBuffer newickString = w.toNewHampshire(p, simpleNewick, writeDistanceToParent);
        return newickString.toString();
    }
    Vector<NJTreeProgressListener> progessListenerVector = new Vector<NJTreeProgressListener>();

    public void addProgessListener(NJTreeProgressListener treeProgessListener) {
        if (treeProgessListener != null) {
            progessListenerVector.add(treeProgessListener);
        }
    }

    public void removeProgessListener(NJTreeProgressListener treeProgessListener) {
        if (treeProgessListener != null) {
            progessListenerVector.remove(treeProgessListener);
        }
    }

    public void broadcastComplete() {
        for (NJTreeProgressListener treeProgressListener : progessListenerVector) {
            treeProgressListener.complete(this);
        }
    }

    public void updateProgress(String state, int percentage) {
        for (NJTreeProgressListener treeProgressListener : progessListenerVector) {
            treeProgressListener.progress(this, state, percentage);
        }
    }

    public void updateProgress(String state, int currentCount, int totalCount) {
        for (NJTreeProgressListener treeProgressListener : progessListenerVector) {
            treeProgressListener.progress(this, state, currentCount, totalCount);
        }
    }

    public static void main(String[] args) {

        try {


            InputStream inStream = TreeConstructor.class.getResourceAsStream("/PF00104_small.fasta");



            FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(inStream, new GenericFastaHeaderParser(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            LinkedHashMap<String,ProteinSequence> proteinSequences = fastaReader.process();
            inStream.close();


            MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound> multipleSequenceAlignment = new MultipleSequenceAlignment<ProteinSequence, AminoAcidCompound>();
            for (ProteinSequence proteinSequence : proteinSequences.values()) {

                multipleSequenceAlignment.addAlignedSequence(proteinSequence);
            }

            long readTime = System.currentTimeMillis();
            TreeConstructor<ProteinSequence, AminoAcidCompound> treeConstructor = new TreeConstructor<ProteinSequence, AminoAcidCompound>(multipleSequenceAlignment, TreeType.NJ, TreeConstructionAlgorithm.PID, new ProgessListenerStub());
            treeConstructor.process();
            long treeTime = System.currentTimeMillis();
            String newick = treeConstructor.getNewickString(true, true);




            System.out.println("Tree time " + (treeTime - readTime));
            System.out.println(newick);

            // treeConstructor.outputPhylipDistances("/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small.fasta.phylip");

        } catch (FileNotFoundException ex) {
            //can't find file specified by args[0]
            ex.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}
