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
 * A TrackLayout implenentation that wraps the sequence smoothly after a set number of residues
 *
 * @author Mark Southern
 * @since 1.5
 */
public class SimpleTrackLayout implements TrackLayout {
    private RangeLocation loc;
    private int wrap = 50;
    private int wrapInc = 10;

    public SimpleTrackLayout() {
    }

    public SimpleTrackLayout(Sequence seq, int wrap) {
        setSequence(seq);
        setWrap(wrap);
    }

    public void setSequence(Sequence seq) {
    }

    public void setRange(RangeLocation loc) {
        this.loc = loc;
    }

    public RangeLocation[] getRanges() {
        int min = loc.getMin();
        int max = loc.getMax();

        int rowCount = ( int ) Math.ceil((( double ) (max - min + 1)) / (( double ) wrap));
        RangeLocation[] ranges = new RangeLocation[ rowCount ];

        int count = 0;

        for (int i = min; i <= max; i += wrap) {
            int newMax = (i + wrap) - 1;

            if (newMax > max) {
                newMax = max;
            }

            RangeLocation loc = new RangeLocation(i, newMax);
            ranges[ count++ ] = loc;
        }

        return ranges;
    }

    public void setWrap(int wrap) {
        this.wrap = wrap;
    }

    public int getWrap() {
        return wrap;
    }

    public int getWrapIncrement() {
        return wrapInc;
    }

    public void setWrapIncrement(int inc) {
        wrapInc = inc;
    }
}
