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

import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;

/** 
 * Encapsulates the basic information you would want from a chromatogram.
 * Read-only.
 *
 *  @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 */
public interface Chromatogram {
    /**
     * The sequence label for the list of called bases.
     */
    public static Object DNA = "dna";

    /**
     * The sequence label for the trace offsets of the called bases.
     */
    public static Object OFFSETS = "trace-offsets";
    
    /** Gets the max intensity from all the traces.  Must be equivalent
     *  to the max of calling {@link #getMax(AtomicSymbol)} on each
     *  of the four non-ambiguous DNA nucleotides.
     *  @return the max intensity 
     */
    public int getMax();
    /** Gets the max intensity on the trace for the specified nucleotide.
     *  @param nucleotide the trace to examine.  Must be a concrete 
     *         (non-ambiguous) nucleotide from the DNA alphabet
     *  @throws IllegalSymbolException when the nucleotide isn't from the DNA
     *          alphabet
     *  @return the max intensity
     */
    public int getMax(AtomicSymbol nucleotide) throws IllegalSymbolException;
    
    /** Returns the length of the trace of the Chromatogram.
     *  @return the number of samples in the trace
     *  @see #getTrace(AtomicSymbol)
     */
    public int getTraceLength();
    /** Returns an array containing the intensities of the sampled waveform
     *  representing the chromatogram trace for base <code>nucleotide</code>.  
     *  This may be a reference the actual internal representation of the 
     *  samples, so callers <b>must not modify it</b>.
     *  <p>
     *  The resulting array for each nucleotide must be {@link #getTraceLength}
     *  <code>int</code>s long.
     *  </p>
     *  @param nucleotide the trace to examine.  Must be the symbol for A, C, G, or T
     *         as provided by {@link org.biojava.bio.seq.DNATools}
     *  @throws IllegalSymbolException if <code>nucleotide</code> isn't in the DNA alphabet
     *  @return an array of integers representing the values of a particular 
     *          chromatogram trace.
     */
    public int[] getTrace(AtomicSymbol nucleotide) throws IllegalSymbolException;
    
    /** 
     * Returns the number of bits of the traces which are significant.  For 
     * instance, if the chromatogram were originally encoded with a single byte 
     * per trace sample, this method must return 8.
     * @return the number of significant bits
     */
    public int getSignificantBits();

    /**
     * Returns an alignment that describes the base calls for this chromatogram.
     * All of the <code>SymbolList</code>s in this alignment must be the same
     * length and that length must equal {@link #getSequenceLength}.
     * <p>
     * The alignment must contain, at the least, two sequences:
     * </p>
     * <ol>
     *   <li>A sequence containing the called bases.  The alphabet of this list
     *       must be {@link org.biojava.bio.seq.DNATools#getDNA()}.  
     *       The label for this list in the alignment must be 
     *       <code>Chromatogram.DNA</code></li>
     *   <li>A sequence containing the trace offsets at which the called bases 
     *       were called.  The alphabet of this list must be an 
     *       {@link org.biojava.bio.symbol.IntegerAlphabet} or a
     *       {@link org.biojava.bio.symbol.IntegerAlphabet.SubIntegerAlphabet}.
     *       The label for this list in the alignment must be 
     *       <code>Chromatogram.OFFSETS</code>.</li>
     * </ol>
     * <p>
     * Implementors may provide other sequences as they see fit.
     * </p>
     * @return an alignment of at least two sequences, as described above.
     */
    public Alignment getBaseCalls();
    /** 
     * Returns the number of bases called by whatever base-calling software
     * analyzed the chromatogram as loaded.  Must equal 
     * <code>{@link #getBaseCalls}.length()</code>.
     * @return the number of bases
     */
    public int getSequenceLength();
    
    /** 
     * Returns a new <code>Chromatogram</code> representing the reverse
     * complement of this one.
     * <p>
     * Implementors should copy the metadata about the chromatogram (i.e., base 
     * calls) as is appropriate to their formats.
     * </p>
     * @return a new chromatogram that is the reverse complement of this one
     */
    public Chromatogram reverseComplement();
}
