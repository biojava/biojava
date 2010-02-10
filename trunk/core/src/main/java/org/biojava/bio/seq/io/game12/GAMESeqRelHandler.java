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
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *  Handles the GAME &lt;<seq_relationship>&gt; element
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMESeqRelHandler
         extends StAXFeatureHandler {
    // <seq_relationship> provide feature positions on
    // specified sequences.
    // This element always has a <span> element that
    // provides the actual coordinates.  This element
    // just adds to it the target sequence name.
    // I will assume there's only one <span> for
    // each <seq_relationship>

    // the only features I will need to consider for
    // now are exons and translate offsets.

    // this is not actually a good place to decide what to 
    // do with incoming data.  I should shove it up the chain
    // to the containing class.

    // conclusion: I think I will forget about exons and
    // and just have all transcripts as compound locations.

    // database columns
    private String type = null;
    private String seq;
//    private String alignment = null;
    Location seqRelLoc;
    StrandedFeature.Strand seqRelStrand;

    // subclass the <span> handler to get access to its data
    private class SpanHandler extends GAMESpanHandler
    {
        private SpanHandler(StAXFeatureHandler staxenv)
        {
            super(staxenv);
//            System.out.println("entering SpanHandler");
        }

        public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler contentHandler) 
        {
            // validate
            super.endElementHandler(nsURI, localName, qName, contentHandler);

            // return the values
//            System.out.println("in SpanHandler: " + loc + " " + strand);
            seqRelLoc = loc;
            seqRelStrand = strand;
        }
    }

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_SEQ_REL_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMESeqRelHandler(staxenv);
            }
        };


    GAMESeqRelHandler(StAXFeatureHandler staxenv) {
        // setup environment
        super(staxenv);

        // set handlers
        // <span> external handler type
        super.addHandler(new ElementRecognizer.ByLocalName("span"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new SpanHandler(staxenv);
                }
            }
        );

        // <alignment>
//        super.addHandler(new ElementRecognizer.ByLocalName("alignment"),
//            new StAXHandlerFactory() {
//                public StAXContentHandler getHandler(StAXFeatureHandler staxenv, long parentID) {
//                    return new AlignmentHandler();
//                }
//            }
//                );
    }


//    private class AlignmentHandler extends StringElementHandlerBase {
//        /**
//         *  Sets the stringValue attribute of the AlignmentHandler object
//         *
//         *@param  s  The new stringValue value
//         */
//        protected void setStringValue(String s) {
//            alignment = s.trim();
//        }
//    }

    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs) {
        // acquire attributes here
        type = attrs.getValue("type");
        seq = attrs.getValue("seq");
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler) 
        throws SAXException 
    {
        // prevalidate
        if ((type == null) || (seq == null) ) {
            System.err.println("improperly formed <span> element.");
        }
    }
}
