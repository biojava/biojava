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

import java.util.Map;

import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.SmallMap;

/**
 * A basic chromatogram implementation which provides public mutators 
 * for setting the various attributes of the chromatogram.  
 * <p>
 * In general, <b>new chromatogram implementations</b> should be derived from
 * {@link AbstractChromatogram}, <b>not this class</b>, as it is generally 
 * undesirable to allow the internal structures of a {@link Chromatogram} to be
 * manipulated externally. This class could still be useful, however, for 
 * programatically generated "chromatograms".
 * </p>
 *
 * @author Rhett Sutphin (<a href="http://genome.uiowa.edu/">UI CBCB</a>)
 * @author Matthew Pocock
 * @since 1.3
 */
public class SimpleChromatogram extends AbstractChromatogram {
    
    /** Creates a new instance of SimpleChromatogram. */
    public SimpleChromatogram() {
        super();
    }
    
    /**
     * Set the DNA and OFFSETS symbol lists for the basecall alignment.  Beware:
     * this method does no consistency checks to be sure that all the offsets
     * are valid indices into the trace arrays.
     *
     * @param dna a symbol list in the DNA alphabet that contains the base calls
     *        for this chromatogram
     * @param offsets a symbol list in an integer or sub-integer alphabet that
     *        contains the locations in the chromatogram for the bases called 
     *        in <code>dna</code>
     * @throws IllegalAlphabetException when the alphabets aren't as specified
     * @throws IllegalArgumentException when the lists aren't the same length
     */
    public void setSymbolLists(SymbolList dna, SymbolList offsets)
    throws IllegalAlphabetException, IllegalArgumentException {
        if (dna.length() != offsets.length())
            throw new IllegalArgumentException("The SymbolLists must be the same length");
        Map map = new SmallMap(2);
        map.put(Chromatogram.DNA, dna);
        map.put(Chromatogram.OFFSETS, offsets);
        super.setBaseCallAlignment(super.createImmutableAlignment(map));
    }
    
    /**
     * Sets the trace array for one of the DNA nucleotides.  The provided
     * array <i>will not be copied</i>, so any modifications to it will be
     * reflected in calls to {@link #getTrace}.
     * <p>
     * If you need to set a new set of traces whose length is different
     * from the old set, you must call {@link #clearTraceValues} first,
     * or you will provoke an <code>IllegalArgumentException</code>.
     * </p>
     *
     * @param nuc the nucleotide for which to set the trace
     * @param trace the sampled intensities along the trace
     * @param maxVal the maximum value on the trace, or {@link java.lang.Integer#MIN_VALUE}
     *        to force this method to calculate it
     * @throws IllegalArgumentException when trace.length is different
     *         from any of the existing (non-null) traces
     * @throws IllegalSymbolException when nuc is not a concrete DNA nucleotide
     */
    public void setTraceValues(AtomicSymbol nuc, int[] trace, int maxVal)
    throws IllegalArgumentException, IllegalSymbolException {
        super.setTrace(nuc, trace, maxVal);
    }
    
    /**
     * Sets all the traces to null.
     */
    public void clearTraceValues() {
        super.clearTraces();
    }
    
    /**
     * Sets the number of significant bits in the data.
     * @param bits how many bits of the trace samples are significant
     * @see Chromatogram#getSignificantBits
     */
    public void setSignificantBits(int bits) {
        super.setBits(bits);
    }
    
    protected AbstractChromatogram reverseComplementInstance() {
        return new SimpleChromatogram();
    }
}
