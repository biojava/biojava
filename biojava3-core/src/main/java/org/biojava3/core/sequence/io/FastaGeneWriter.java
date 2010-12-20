/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 * Created on 01-21-2010
 */

package org.biojava3.core.sequence.io;

import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;

import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ChromosomeSequence;
import org.biojava3.core.sequence.ExonSequence;
import org.biojava3.core.sequence.GeneSequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.template.FastaHeaderFormatInterface;

/**
 * A Gene sequence has a Positive or Negative Strand where we want to write out to a stream the 5 to 3 prime version.
 * It is also an option to write out the gene sequence where the exon regions are upper case
 * 6/22/2010 FastaWriter needs to be sequence aware to handle writing out a GeneSequence which is negative Strand with the proper sequence
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaGeneWriter {

    boolean showExonUppercase = false;
    OutputStream os;
    Collection<GeneSequence> sequences;
    FastaHeaderFormatInterface<GeneSequence, NucleotideCompound> headerFormat;
    private int lineLength = 60;
/**
 *
 * @param os
 * @param sequences
 * @param headerFormat
 * @param showExonUppercase
 */
    public FastaGeneWriter(OutputStream os, Collection<GeneSequence> sequences, FastaHeaderFormatInterface<GeneSequence, NucleotideCompound> headerFormat, boolean showExonUppercase) {
        this(os, sequences, headerFormat, showExonUppercase, 60);
    }
/**
 *
 * @param os
 * @param sequences
 * @param headerFormat
 * @param showExonUppercase
 * @param lineLength
 */
    public FastaGeneWriter(OutputStream os, Collection<GeneSequence> sequences, FastaHeaderFormatInterface<GeneSequence, NucleotideCompound> headerFormat, boolean showExonUppercase, int lineLength) {
        this.os = os;
        this.sequences = sequences;
        this.headerFormat = headerFormat;
        this.lineLength = lineLength;
        this.showExonUppercase = showExonUppercase;
    }
/**
 *
 * @throws Exception
 */
    public void process() throws Exception {
        byte[] lineSep = System.getProperty("line.separator").getBytes();

        for (GeneSequence sequence : sequences) {
            String header = headerFormat.getHeader(sequence);
            os.write('>');
            os.write(header.getBytes());
            os.write(lineSep);

            int compoundCount = 0;
            String seq = "";
            //GeneSequence currently has a strand attribute to indicate direction

            seq = sequence.getSequence5PrimeTo3Prime().getSequenceAsString();
            if (showExonUppercase) {
                StringBuilder sb = new StringBuilder(seq.toLowerCase());
                int geneBioBegin = sequence.getBioBegin();
                int geneBioEnd = sequence.getBioEnd();
                for (ExonSequence exonSequence : sequence.getExonSequences()) {
                    int featureBioBegin = 0;
                    int featureBioEnd = 0;
                    if (sequence.getStrand() != Strand.NEGATIVE) {
                        featureBioBegin = exonSequence.getBioBegin() - geneBioBegin;
                        featureBioEnd = exonSequence.getBioEnd() - geneBioBegin;
                    } else {
                        featureBioBegin = geneBioEnd - exonSequence.getBioEnd();
                        featureBioEnd = geneBioEnd - exonSequence.getBioBegin();
                    }
                    if (featureBioBegin < 0 || featureBioEnd < 0 || featureBioEnd > sb.length() || featureBioBegin > sb.length()) {
                        System.out.println("Bad Feature " + sequence.getAccession().toString() + " " + sequence.getStrand() + " " + geneBioBegin + " " + geneBioEnd + " " + exonSequence.getBioBegin() + " " + exonSequence.getBioEnd());
                    } else {
                        for (int i = featureBioBegin; i <= featureBioEnd; i++) {
                            char ch = sb.charAt(i);
                            //probably not the fastest but the safest way if language is not standard ASCII
                            String temp = ch + "";
                            ch = temp.toUpperCase().charAt(0);
                            sb.setCharAt(i, ch);
                        }
                    }
                }
                seq = sb.toString();
            }

            for (int i = 0; i < seq.length(); i++) {
                os.write(seq.charAt(i));
                compoundCount++;
                if (compoundCount == lineLength) {
                    os.write(lineSep);
                    compoundCount = 0;
                }

            }


            //If we had sequence which was a reciprocal of line length
            //then don't write the line terminator as this has already written
            //it
            if ((sequence.getLength() % getLineLength()) != 0) {
                os.write(lineSep);
            }
        }
    }

    /**
     * @return the lineLength
     */
    public int getLineLength() {
        return lineLength;
    }

    /**
     * @param lineLength the lineLength to set
     */
    public void setLineLength(int lineLength) {
        this.lineLength = lineLength;
    }

    public static void main(String[] args) {

        try {
            ArrayList<GeneSequence> sequences = new ArrayList<GeneSequence>();
            ChromosomeSequence seq1 = new ChromosomeSequence("ATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGC");
            GeneSequence gene1 = seq1.addGene(new AccessionID("gene1"), 1, 20, Strand.POSITIVE);

            gene1.addExon(new AccessionID("t1_1_10"), 1, 10);
            gene1.addExon(new AccessionID("t1_12_15"), 12, 15);
            GeneSequence gene2 = seq1.addGene(new AccessionID("gene2"), 1, 20, Strand.NEGATIVE);

            gene2.addExon(new AccessionID("t2_1_10"), 1, 10);
            gene2.addExon(new AccessionID("t2_12_15"), 12, 15);
            sequences.add(gene1);
            sequences.add(gene2);


            FastaGeneWriter fastaWriter = new FastaGeneWriter(System.out, sequences, new GenericFastaHeaderFormat<GeneSequence, NucleotideCompound>(), true);
            fastaWriter.process();


        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
