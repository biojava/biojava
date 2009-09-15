/**
 *  BioJava development code This code may be freely distributed and modified
 *  under the terms of the GNU Lesser General Public Licence. This should be
 *  distributed with the code. If you do not have a copy, see:
 *  http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 *  jointly by the individual authors. These should be listed in
 *
 *@author    doc comments. For more information on the BioJava project and its
 *      aims, or to join the biojava-l mailing list, visit the home page at:
 *      http://www.biojava.org/
 */

package org.biojava.bio.seq.io.game12;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.stax.IntElementHandlerBase;
import org.biojava.utils.stax.StAXContentHandler;

/**
 *  Handles the GAME &lt;<span>&gt; element.
 *  Subclass this to parse &lt;<span>&gt; and get your result somewhere useful.
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMESpanHandler
         extends StAXFeatureHandler {
    // <span> captures sequence locations.

    // database columns
    int start;
    int end;
    boolean gotStart = false;
    boolean gotEnd = false;
    Location loc;
    StrandedFeature.Strand strand;

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_SPAN_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMESpanHandler(staxenv);
            }
        };


    /**
     *  Constructor for the GAMESpanHandler object
     *
     *@param  staxenv   Description of the Parameter
     *@param  parentID  Description of the Parameter
     */
    GAMESpanHandler(StAXFeatureHandler staxenv) {
        // setup environment
        super(staxenv);

        // setup handlers
        // <start>
        super.addHandler(new ElementRecognizer.ByLocalName("start"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new StartHandler();
                }
            }
                );
        // <end>
        super.addHandler(new ElementRecognizer.ByLocalName("end"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new EndHandler();
                }
            }
                );
    }

    private class StartHandler extends IntElementHandlerBase {
        /**
         *  Sets the intValue attribute of the StartHandler object
         *
         *@param  startVal  The new intValue value
         */
        protected void setIntValue(int startVal) {
            start = startVal;
            gotStart = true;
        }
    }


    private class EndHandler extends IntElementHandlerBase {
        /**
         *  Sets the intValue attribute of the EndHandler object
         *
         *@param  endVal  The new intValue value
         */
        protected void setIntValue(int endVal) {
            end = endVal;
            gotEnd = true;
        }
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler) {
        // prevalidate
        if (!gotStart || !gotEnd) {
            return;
        }

        // create a RangeLocation that embodies the info
        // remember that in their nomenclature, a point
        // location with strand can only be done by
        // two coordinates separated by one.

        if (start < end) strand = StrandedFeature.POSITIVE;
        else if (start > end) strand = StrandedFeature.NEGATIVE;
        else strand = StrandedFeature.UNKNOWN;

        int min = Math.min(start,end);
        int max = Math.max(start,end);

        if (Math.abs(start - end) == 1) {
            // they've got a point location in mind
            loc = new PointLocation(min);    
        }
        else {
            // range location required
            loc = new RangeLocation(min, max);
        }
    }
}

