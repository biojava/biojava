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
 */

package org.biojava.bio.chromatogram;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;

import org.biojava.bio.BioError;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolListViews;
import org.biojava.utils.ChangeListener;

/**
 * A basic, abstract implementation of {@link Chromatogram}.  Provides
 * protected setters so that subclasses may set the value of the various
 * properties of a chromatogram.
 *
 * Chromatograms should be created using {@link ChromatogramFactory} or a
 * parser for a particular file format.
 * 
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 *
 * 
 */
public abstract class AbstractChromatogram implements Chromatogram {
    private static int A = 0; // Array index for 'A'-related material
    private static int C = 1; // Array index for 'C'-related material
    private static int G = 2; // Array index for 'G'-related material
    private static int T = 3; // Array index for 'T'-related material

    /** 2D array containing all the samples for all 4 traces. */
    private int[][] traceSample;
    /** Array containing the single highest value for each trace. */
    private int[] maxTraceValue;
    /** The immutable Alignment that will be returned by getBaseCalls(). */
    private Alignment baseCalls;
    private int significantBits;

  /**
   * Create a new AbstractChromatogram.
   */
    public AbstractChromatogram() {
        maxTraceValue = new int[4];
        traceSample = new int[4][];
        significantBits = 0;
    }

    public int[] getTrace(AtomicSymbol nucleotide)
    throws IllegalSymbolException {
        return traceSample[nucToIndex(nucleotide)];
    }

    public int getTraceLength() {
        return traceSample[0].length;
    }

    public int getMax() {
        try {
            return Math.max(
                    Math.max(getMax(DNATools.a()), getMax(DNATools.c())),
                    Math.max(getMax(DNATools.g()), getMax(DNATools.t()))
                   );
        } catch (IllegalSymbolException ise) {
            throw new BioError("Can't happen", ise);
        }
    }

    public int getMax(AtomicSymbol nucleotide) throws IllegalSymbolException {
        return maxTraceValue[nucToIndex(nucleotide)];
    }

  /**
   * Return the total number of base calls.
   *
   * @return the total number of base calls
   */
    public Alignment getBaseCalls()  { return baseCalls; }

  /**
   * Return the sequence length.
   *
   * @return the sequence length
   */
    public int getSequenceLength()   { return baseCalls.length(); }

  /**
   * Return the number of significant bits.
   *
   * @return the significant bits
   */
    public int getSignificantBits()  { return significantBits; }

    /**
     * Convert a DNA symbol to an internal array index.
     *
     * @param nuc  the nucleotide to convert
     * @return an integer giving it's internal index value
     * @throws IllegalSymbolException if the symbol is not recognised
     */
    private final int nucToIndex(AtomicSymbol nuc) throws IllegalSymbolException {
        if      (nuc == DNATools.a()) return A;
        else if (nuc == DNATools.c()) return C;
        else if (nuc == DNATools.g()) return G;
        else if (nuc == DNATools.t()) return T;
        else
            throw new IllegalSymbolException("The symbol " +
                nuc.getName() + " (" +
                nuc.getClass().getName() +
                ") is not in the DNA alphabet");
    }

    /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
     * Protected mutators for subclasses to initialize fields
     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

    /**
     * Provides the list of base calls.
     * @param align the base call alignment
     * @throws NoSuchElementException when align doesn't contain alignments with
     *         the required DNA and OFFSETS labels
     * @throws IllegalArgumentException the lists in align aren't all the same length.
     * @throws IllegalAlphabetException if the required lists don't have the
     *         correct alphabets.  See the documentation of
     *         {@link Chromatogram#getBaseCalls() } for details.
     * @see Chromatogram#getBaseCalls
     */
    protected final void setBaseCallAlignment(Alignment align)
    throws IllegalAlphabetException, IllegalArgumentException,
           NoSuchElementException {
        SymbolList dna, offsets;
        try {
            dna = align.symbolListForLabel(Chromatogram.DNA);
            offsets = align.symbolListForLabel(Chromatogram.OFFSETS);
        } catch (NoSuchElementException nsee) {
            throw nsee;
        }
        if (dna.getAlphabet() != DNATools.getDNA())
            throw new IllegalAlphabetException("DNA list has inappropriate alphabet");
        if (! (offsets.getAlphabet() instanceof IntegerAlphabet ||
               offsets.getAlphabet() instanceof IntegerAlphabet.SubIntegerAlphabet))
            throw new IllegalAlphabetException("Offsets list has inappropriate alphabet");

        baseCalls = align;
    }

    /**
     * Sets the trace data structures to null.  If a subclass needs to replace
     * the traces with new traces of a different length, this method must be
     * called first to avoid provoking an IllegalArgumentException from setTrace
     * due to a length mismatch.
     */
    protected final void clearTraces() {
        maxTraceValue[A] = maxTraceValue[C] = maxTraceValue[G]
                         = maxTraceValue[T] = Integer.MIN_VALUE;
        traceSample[A] = traceSample[C] = traceSample[G]
                       = traceSample[T] = null;
    }

    /**
     * Provides the trace samples for a particular nucleotide.
     * @param nuc A DNA nucleotide
     * @param trace the trace samples themselves
     * @param maxVal the maximum value in the trace array.  If this value
     *        is {@link Integer#MIN_VALUE}, this method will do a linear
     *        search of trace to determine the max.
     * @throws IllegalArgumentException when trace.length is different
     *         from any of the existing (non-null) traces
     * @throws IllegalSymbolException when nuc is not a concrete DNA nucleotide
     */
    protected final void setTrace(AtomicSymbol nuc, int[] trace, int maxVal)
    throws IllegalArgumentException, IllegalSymbolException {
        int idx = nucToIndex(nuc);
        if (trace == null) {
            traceSample[idx] = null;
            maxTraceValue[idx] = Integer.MIN_VALUE;
        }
        else {
            // check for length mismatches
            for (int i = 0 ; i < 4 ; i++) {
                if (traceSample[i] != null && traceSample[i].length != trace.length) {
                    throw new IllegalArgumentException(
                        "All traces must be the same length.  " + trace.length  +
                        " != " + traceSample[i].length);
                }
            }
            // special case to do linear search for max
            if (maxVal == Integer.MIN_VALUE) {
                maxVal = trace[0];
                for (int i = 1 ; i < trace.length ; i++)
                    maxVal = Math.max(maxVal, trace[i]);
            }
            traceSample[idx] = trace;
            maxTraceValue[idx] = maxVal;
        }
    }

    /** Sets the number of significant bits in the trace samples.
     *  @param bits a non-negative integer indicating the number of
     *         significant bits in each trace sample
     *  @throws IllegalArgumentException when <code>bits</code> is negative
     */
    protected final void setBits(int bits)
    throws IllegalArgumentException {
        if (bits < 0)
            throw new IllegalArgumentException("Invalid number of significant bits");
        this.significantBits = bits;
    }

    /** Returns a new instance of this AbstractChromatogram subclass for use in
     * {@link #reverseComplement}.
     *
     * @return a reverse-complemented AbstractChromatogram
     */
    protected abstract AbstractChromatogram reverseComplementInstance();

    public Chromatogram reverseComplement() {
        AbstractChromatogram rev = this.reverseComplementInstance();
        try {
            rev.setTrace(DNATools.a(), reverse(traceSample[T]), maxTraceValue[T]);
            rev.setTrace(DNATools.c(), reverse(traceSample[G]), maxTraceValue[G]);
            rev.setTrace(DNATools.g(), reverse(traceSample[C]), maxTraceValue[C]);
            rev.setTrace(DNATools.t(), reverse(traceSample[A]), maxTraceValue[A]);
        } catch (IllegalSymbolException ise) {
            throw new BioError( "Can't happen -- all symbols are explicit and legal", ise);
        }
        Alignment revBC = reverseComplementBaseCalls();
        try {
            rev.setBaseCallAlignment(revBC);
        } catch (IllegalAlphabetException iae) {
            throw new BioError("Can't happen unless reverseComplementBaseCalls or reverseComplementBaseCallList have been overridden out-of-spec", iae);
        } catch (NoSuchElementException nsee) {
            throw new BioError("Can't happen unless reverseComplementBaseCalls or reverseComplementBaseCallList have been overridden out-of-spec", nsee);
        } catch (IllegalArgumentException iae) {
            throw new BioError("Can't happen unless reverseComplementBaseCalls or reverseComplementBaseCallList have been overridden out-of-spec", iae);
        }
        rev.setBits(this.significantBits);
        return rev;
    }

    /**
     * Returns a new base call alignment that is the reverse complement of
     * one in this chromatogram.  This is achieved by calling
     * {@link #reverseComplementBaseCallList} for each label in the current
     * base call alignment.  When that method returns null, no list will
     * appear in reverse complement base call alignment with the null-provoking
     * label.  For this reason, subclasses are encouraged to override
     * {@link #reverseComplementBaseCallList} to handle any additional per-base
     * metadata that they store.
     * <p>
     * This implementation should be safely inheritable
     * for all chromatogram implementations, unless one just doesn't want
     * base calls on its reverse complement output.  If this is the case,
     * it should override this method to return null.
     * </p>
     *
     * @return a new {@link Alignment} that is the reverse complement of the
     *          one in the current chromatogram
     */
    protected Alignment reverseComplementBaseCalls() {
        Alignment current = this.getBaseCalls();
        Map revBCMap = new HashMap();
        Object curLabel;
        SymbolList revList;
        for (Iterator it = current.getLabels().iterator() ; it.hasNext() ;) {
            curLabel = it.next();
            revList = reverseComplementBaseCallList(curLabel);
            if (revList != null)
                revBCMap.put(curLabel, revList);
        }
        return createImmutableAlignment(revBCMap);
    }

    /**
     * Return a symbol list containing the reverse complement of the base call
     * data for the given label.  The returned list will be stored in the
     * reverse complement's base call alignment under the same label.
     * <p>
     * Implementation note: subclasses which do not use an {@link IntegerAlphabet} for
     * their offsets lists must override this method, at least for the case where
     * <code>label == {@link #OFFSETS}</code>.
     * </p>
     *
     * @param label the label Object
     * @return an appropriately reverse-complemented SymbolList, or null if the
     *          label is unhandled.
     */
    protected SymbolList reverseComplementBaseCallList(Object label) {
        if (label == DNA) {
            try {
                return SymbolListViews.reverse(DNATools.complement(this.getBaseCalls().symbolListForLabel(DNA)));
            } catch (IllegalAlphabetException iae) {
                throw new BioError("Can't happen unless the DNA list has been set improperly", iae);
            }
        }
        else if (label == OFFSETS) {
            SymbolList curOffsets = this.getBaseCalls().symbolListForLabel(OFFSETS);
            List revOffsets = new ArrayList(curOffsets.length());
            IntegerAlphabet alpha = (IntegerAlphabet) curOffsets.getAlphabet();
            IntegerAlphabet.IntegerSymbol sym;
            for (int i = curOffsets.length() ; i > 0 ; i--) {
                sym = (IntegerAlphabet.IntegerSymbol) curOffsets.symbolAt(i);
                revOffsets.add(alpha.getSymbol(this.getTraceLength() - sym.intValue()));
            }
            try {
                return createImmutableSymbolList(alpha, revOffsets);
            } catch (IllegalSymbolException ise) {
                throw new BioError("Can't happen -- revOffsets was just created with only symbols from alpha");
            }
        }
        else {
            return null;
        }
    }

    /**
     * A factory method for creating new symbol lists with a given alphabet.
     * The default implementation should be fine for nearly all cases, but the
     * option is given in case a developer wishes to use a more memory
     * efficient implementation.
     *
     * @param alpha the {@link Alphabet} for the new list
     * @param syms the symbols to put in the new list
     * @throws IllegalSymbolException when alpha and syms are incompatible
     * @throws ClassCastException when any object in syms isn't a Symbol
     * @return a new immutable SymbolList containing all the given symbols using
     *          the given Alphabet
     */
    protected SymbolList createImmutableSymbolList(Alphabet alpha, List syms)
    throws IllegalSymbolException, ClassCastException {
        SymbolList symlist = new SimpleSymbolList(alpha, syms);
        symlist.addChangeListener(ChangeListener.ALWAYS_VETO);
        return symlist;
    }

    /**
     * A factory method for creating new immutable alignments, particularly
     * for use as base call alignments.  The default implementation should
     * be fine for nearly all cases.
     *
     * @param labelsToSymLists a {@link Map} whose keys are desired labels
     *        for the alignment and whose values are the SymbolLists.
     *        All the SymbolLists must be the same length.
     * @return a new Alignment
     * @throws IllegalArgumentException if the lists aren't all the same length
     * @throws ClassCastException if any of the values in the map aren't
     *         SymbolLists
     */
    protected Alignment createImmutableAlignment(Map labelsToSymLists)
    throws IllegalArgumentException, ClassCastException {
        Alignment sa = new SimpleAlignment(labelsToSymLists);
        sa.addChangeListener(ChangeListener.ALWAYS_VETO);
        return sa;
    }

    /**
     *  Utility method for reversing an int[] array. Visible for subclass use.
     *
     * @param src  the source array
     * @return an array of the same length containing all values in src in
     *     reverse order
     */
    protected final static int[] reverse(int[] src) {
        int[] dst = new int[src.length];
        for (int i = 0 ; i < src.length ; i++)
            dst[src.length - i - 1] = src[i];
        return dst;
    }
}
