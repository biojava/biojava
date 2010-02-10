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

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;


/**
 * An implementation of TrackLayout that that wraps a sequence over an arbitrary set of ranges
 *
 * @author Mark Southern
 * @since 1.5
 */
public class UserDefinedTrackLayout implements TrackLayout {
    private RangeLocation globalLocation;
    private RangeLocation[] ranges;
    private int wrap = 6;
    private int wrapInc = 1;

    public UserDefinedTrackLayout(RangeLocation[] ranges) {
        this.ranges = ranges;
    }

    public void setSequence(Sequence seq) {
    }

    public void setRange(RangeLocation loc) {
        this.globalLocation = loc;
    }

    public void setWrap(int wrap) {
        this.wrap = wrap;
    }

    public int getWrap() {
        return wrap;
    }

    public RangeLocation[] getRanges() {
        List list = new ArrayList();

        for (int i = 0; i < ranges.length; i++) {
            Location loc = LocationTools.intersection(globalLocation, ranges[ i ]);

            if (! loc.equals(Location.empty)) {
                list.add(loc);
            }
        }

        return ( RangeLocation[] ) list.toArray(new RangeLocation[ 0 ]);
    }

    public int getWrapIncrement() {
        return wrapInc;
    }

    public void setWrapIncrement(int inc) {
        wrapInc = inc;
    }
}
