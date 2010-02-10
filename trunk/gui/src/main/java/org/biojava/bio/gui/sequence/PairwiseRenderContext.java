/*
 * BioJava development code
 * 
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 * 
 * http://www.gnu.org/copyleft/lesser.html
 * 
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 * 
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 * 
 * http://www.biojava.org
 * 
 */

package org.biojava.bio.gui.sequence;

import java.awt.Point;

import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * <code>PairwiseRenderContext</code> encapsulates information
 * required for the rendering of a pair of sequences. No assumption is
 * made as to whether the sequences are to be rendered in different
 * directions (as in a dotplot) or in the same direction; this is left
 * to the implementation. The leading and trailing borders refer to
 * the primary sequence only.
 *
 * @author Keith James
 * @since 1.2
 */
public interface PairwiseRenderContext extends SequenceRenderContext
{
    /**
     * <code>getSecondaryDirection</code> returns the direction in
     * which the secondary sequence is rendered. This may be either
     * HORIZONTAL or VERTICAL.
     *
     * @return an <code>int</code>.
     */
    public int getSecondaryDirection();

    /**
     * <code>getSecondarySymbols</code> returns the symbols of the
     * secondary sequence.
     *
     * @return a <code>SymbolList</code>.
     */
    public SymbolList getSecondarySymbols();

    /**
     * <code>getSecondaryFeatures</code> returns the features on the
     * secondary sequence.
     *
     * @return a <code>FeatureHolder</code>.
     */
    public FeatureHolder getSecondaryFeatures();

    /**
     * <code>getSecondaryRange</code> returns the range of the
     * secondary sequence currently rendered.
     *
     * @return a <code>RangeLocation</code>.
     */
    public RangeLocation getSecondaryRange();

    /**
     * <code>secondarySequenceToGraphics</code> converts a sequence
     * coordinate on the secondary sequence to a graphical position.
     *
     * @param sequencePos an <code>int</code>.
     *
     * @return a <code>double</code>.
     */
    public double secondarySequenceToGraphics(int sequencePos);

    /**
     * <code>graphicsToSecondarySequence</code> converts a graphical
     * position to a sequence coordinate on the secondary sequence.
     *
     * @param graphicsPos a <code>double</code>.
     *
     * @return an <code>int</code>.
     */
    public int graphicsToSecondarySequence(double graphicsPos);

    /**
     * <code>graphicsToSecondarySequence</code> converts a graphical
     * position to a secondary sequence index.
     *
     * @param point a <code>Point</code>.
     *
     * @return an <code>int</code>.
     */
    public int graphicsToSecondarySequence(Point point);
}
