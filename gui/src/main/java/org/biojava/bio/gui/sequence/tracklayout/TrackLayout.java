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
package org.biojava.bio.gui.sequence.tracklayout;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.RangeLocation;


/**
 * An interface for the handling of the layout of a WrappedSequencePanel.
 *
 * @author Mark Southern
 * @since 1.5
 */
public interface TrackLayout {
    void setSequence(Sequence seq);

    void setRange(RangeLocation loc);

    RangeLocation[] getRanges();

    void setWrap(int wrap);

    int getWrap();

    void setWrapIncrement(int inc);

    int getWrapIncrement();
}
