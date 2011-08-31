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
 * Created on DATE
 *
 */

package org.biojava3.core.sequence;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.biojava3.core.sequence.io.util.IOUtils;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.LightweightProfile;
import org.biojava3.core.sequence.template.Sequence;

/**
 * Implements a minimal data structure for reading and writing a sequence alignment.  The full {@code Profile} data
 * structure in the alignment module provides additional functionality.
 *
 * @author Scooter Willis
 * @author Mark Chapman
 */
public class MultipleSequenceAlignment<S extends Sequence<C>, C extends Compound> implements LightweightProfile<S, C> {

    private List<S> sequences = new ArrayList<S>();
    private Integer length = null;

    /**
     * A sequence that has been aligned to other sequences will have inserts. 
     * @param sequence
     */
    public void addAlignedSequence(S sequence){
        if(length == null){
            length = sequence.getLength();
        }
        if(sequence.getLength() != length){
            throw new IllegalArgumentException(sequence.getAccession() + " length = " + sequence.getLength() +
                    " not equal to MSA length = " + length);
        }
        sequences.add(sequence);
    }

    /**
     * Remove a sequence
     * @param sequence
     * @return flag
     */
    public boolean removeAlignedSequence(S sequence){
        return sequences.remove(sequence);
    }
//methods for LightweightProfile

    /**
     * Uses bioIndex starting at 1 instead of 0
     * @param listIndex
     * @return sequence
     */
     

    @Override
    public S getAlignedSequence(int listIndex) {
        return sequences.get(listIndex - 1);
    }

    /**
     * Get the list of sequences
     * @return list of sequences
     */
    @Override
    public List<S> getAlignedSequences() {
        return Collections.unmodifiableList(sequences);
    }

    /**
     * Get a list of compounds at a sequence position
     * @param alignmentIndex
     * @return compounds
     */
    @Override
    public List<C> getCompoundsAt(int alignmentIndex) {
        List<C> column = new ArrayList<C>();
        for (S s : sequences) {
            column.add(s.getCompoundAt(alignmentIndex));
        }
        return Collections.unmodifiableList(column);
    }

    /**
     * Get the Compounds defined in the first sequence
     * @return get compound set
     */
    @Override
    public CompoundSet<C> getCompoundSet() {
        return sequences.get(0).getCompoundSet();
    }

    /**
     * Get the length of the MSA where it is assumed that
     * all sequence position
     * @return length of MSA
     */
    @Override
    public int getLength() {
        return length;
    }

    /**
     * Get the number of sequences in the MSA
     * @return nr of sequences
     */
    @Override
    public int getSize() {
        return sequences.size();
    }

    /**
     * Get a string representation of the MSA with a fixed width
     * @param width
     * @return String
     */
    @Override
    public String toString(int width) {
        return toString(width, null, IOUtils.getIDFormat(sequences), true, true, true, false);
    }

    /**
     * Support for different MSA formats
     * @param format
     * @return String in one of the supported file formats.
     */
    @Override
    public String toString(StringFormat format) {
        switch (format) {
        case ALN:
        case CLUSTALW:
        default:
            return toString(60, String.format("CLUSTAL W MSA from BioJava%n%n"), IOUtils.getIDFormat(sequences) +
                    "   ", true, false, true, false);
        case FASTA:
            return toString(60, null, ">%s%n", false, false, false, false);
        case GCG:
        case MSF:
            return toString(50, IOUtils.getGCGHeader(sequences), IOUtils.getIDFormat(sequences), true, false, false,
                    false);
        case PDBWEB:
            return toString(60, null, "%s", true, false, true, true);
        }
    }

    /**
     * String representation of the MSA
     * @return String
     */

    @Override
    public String toString() {
        return toString(getLength(), null, null, false, false, false, false);
    }

    // helper methods

    /**
     * Helper method that does all the formating work
     * @param width
     * @param header
     * @param idFormat
     * @param interlaced
     * @param aligIndices
     * @param aligConservation
     * @param webDisplay
     * @return String
     */
    // creates formatted String
    private String toString(int width, String header, String idFormat, boolean interlaced, boolean aligIndices,
            boolean aligConservation, boolean webDisplay) {

        // TODO handle circular alignments
        StringBuilder s = (header == null) ? new StringBuilder() : new StringBuilder(header);

        if (webDisplay && sequences.size() == 2) {
            s.append("<div><pre>");
        }

        width = Math.max(1, width);
        if (interlaced) {
            String aligIndFormat = "%-" + Math.max(1, width / 2) + "d %" + Math.max(1, width - (width / 2) - 1) +
                    "d%n";
            for (int i = 0; i < getLength(); i += width) {
                int start = i + 1, end = Math.min(getLength(), i + width);
                if (i > 0) {
                    s.append(String.format("%n"));
                }
                if (aligIndices) {
                    if (end < i + width) {
                        int line = end - start + 1;
                        aligIndFormat = "%-" + Math.max(1, line / 2) + "d %" + Math.max(1, line - (line / 2) - 1) +
                                "d%n";
                    }
                    if (idFormat != null) {
                        s.append(String.format(idFormat, ""));
                    }
                    s.append(String.format(aligIndFormat, start, end));
                }
                int counter = 0;
                for (S as : sequences) {
                    counter++;
                    if (webDisplay && sequences.size() == 2) {
                        printSequenceAlignmentWeb(s, counter, idFormat, start, end);
                    } else {
                        if (idFormat != null) {
                            s.append(String.format(idFormat, as.getAccession()));
                        }
                        s.append(as.getSubSequence(start, end).getSequenceAsString());
                        s.append(String.format("%n"));
                    }
                    if (aligConservation && sequences.size() == 2 && counter == 1) {
                        printConservation(s, idFormat, start, end, webDisplay);
                    }
                }
            }
        } else {
            for (S as : sequences) {
                if (idFormat != null) {
                    s.append(String.format(idFormat, as.getAccession()));
                }
                for (int i = 0; i < getLength(); i += width) {
                    int start = i + 1, end = Math.min(getLength(), i + width);
                    s.append(as.getSubSequence(start, end).getSequenceAsString());
                    s.append(String.format("%n"));
                }
            }
        }

        if (webDisplay && aligConservation && sequences.size() == 2) {
            s.append(IOUtils.getPDBLegend());
        }
        return s.toString();
    }

    /**
     *
     * @param s
     * @param counter
     * @param idFormat
     * @param start
     * @param end
     */
    private void printSequenceAlignmentWeb(StringBuilder s, int counter, String idFormat, int start, int end) {
        S as = sequences.get(counter - 1), seq1 = sequences.get(0), seq2 = sequences.get(1);

        if (idFormat != null) {
            s.append(String.format(idFormat, as.getAccession()));
        }

        String mySeq = as.getSubSequence(start, end).getSequenceAsString();
        String s1 = seq1.getSubSequence(start, end).getSequenceAsString();
        String s2 = seq2.getSubSequence(start, end).getSequenceAsString();
        CompoundSet<C> cs = getCompoundSet();

        for (int i = 0; i < s1.length(); i++) {
            if (i >= s2.length() || i >= mySeq.length())
                break;
            char c1 = s1.charAt(i);
            char c2 = s2.charAt(i);
            char c = mySeq.charAt(i);
            s.append(IOUtils.getPDBCharacter(true, c1, c2, cs.compoundsEquivalent(seq1.getCompoundAt(i),
                    seq2.getCompoundAt(i)), c));
        }

        s.append(String.format("%n"));
    }

    /**
     *
     * @param s
     * @param idFormat
     * @param start
     * @param end
     * @param webDisplay
     */
    private void printConservation(StringBuilder s, String idFormat, int start, int end, boolean webDisplay) {
        S seq1 = sequences.get(0), seq2 = sequences.get(1);

        if (idFormat != null) {
            AccessionID ac1 = sequences.get(0).getAccession();
            String id1 = (ac1 == null) ? "null" : ac1.getID();
            id1 = id1.replaceAll(".", " ");
            s.append(String.format(idFormat, id1));
        }

        String s1 = seq1.getSubSequence(start, end).getSequenceAsString();
        String s2 = seq2.getSubSequence(start, end).getSequenceAsString();
        CompoundSet<C> cs = getCompoundSet();

        for (int i = 0; i < s1.length(); i++) {
            if (i >= s2.length())
                break;
            char c1 = s1.charAt(i);
            char c2 = s2.charAt(i);
            s.append(IOUtils.getPDBConservation(webDisplay, c1, c2, cs.compoundsEquivalent(seq1.getCompoundAt(i),
                    seq2.getCompoundAt(i))));
        }

        s.append(String.format("%n"));
    }

}
