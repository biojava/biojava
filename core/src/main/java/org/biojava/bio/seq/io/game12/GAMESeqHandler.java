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

import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *  Handles the GAME &lt;seq&gt; element
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMESeqHandler
         extends StAXFeatureHandler {
    // the <seq> element supplies clone name and length.
    // other data includes a description of the sequence.
    // <seq> does not necessarily have an accompanying <residues>
    // seq appears to be near terminal in that it doesn't contain
    // complex structures with their own unique_ids within it.
    // this means that I can safely omit the duplicate entry without
    // risking that some nested element within it goes wrong.

    //database columns
    private String seqId;

    private int seqLength = 0;
//    private byte[] buffer = null;
//    private boolean hasResidues = false;
//    private boolean nonUniqueEntry = false;

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_SEQ_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMESeqHandler(staxenv);
            }
        };


    /**
     *  Constructor for the GAMESeqHandler object
     *
     *@param  staxenv   Description of the Parameter
     *@param  parentID  Description of the Parameter
     */
    GAMESeqHandler(StAXFeatureHandler staxenv) {
        // stash environment
        super(staxenv);

        // setup handlers
        // <name>
        super.addHandler(new ElementRecognizer.ByLocalName("name"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new NameHandler();
                }
            }
                );

        // <description>
        super.addHandler(new ElementRecognizer.ByLocalName("description"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new DescriptionHandler();
                }
            }
                );
        // <residues>
//        super.addHandler(new ElementRecognizer.ByLocalName("residues"),
//                GAMEResiduesHandler.GAME_RESIDUES_HANDLER_FACTORY);
    }


    /**
     *  Description of the Class
     *
     *@author     david
     *@created    19 January 2002
     */
    private class NameHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the NameHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
        }
    }


    /**
     *  Description of the Class
     *
     *@author     david
     *@created    19 January 2002
     */
    private class DescriptionHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the DescriptionHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
        }
    }


    /**
     *  Gets the sequenceLength attribute of the GAMESeqHandler object
     *
     *@return    The sequenceLength value
     */
    public int getSequenceLength() {
        return seqLength;
    }


//    /**
//     *  Description of the Method
//     *
//     *@param  buffer  Description of the Parameter
//     */
//    public void returnSequenceBuffer(byte[] buffer) {
//        this.buffer = buffer;
//    }


    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs)
             throws SAXException 
    {
        // pick up attributes
        seqId = attrs.getValue("id").trim();        
        String length = attrs.getValue("length").trim();        

        // return results
        try {
            listener.setName(seqId);
            listener.addSequenceProperty("length", length);
        }
        catch (ParseException pe) {
            throw new SAXException("could not set sequence properties.");
        }
    }


    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler staxHandler)
             throws SAXException {
    }
}
