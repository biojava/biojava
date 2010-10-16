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
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *  Handles the root GAME element
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMEHandler
         extends StAXFeatureHandler {
    // there is only one GAME element encompassing the entire file
    // in Gadfly GAME files.

    /**
     *  Constructor for the GAMEHandler object
     *
     */
    public GAMEHandler(SeqIOListener listener) {
        // initialise environment
        this.staxenv = this;
        this.listener = listener;

        // setup handlers
        // <seq>
        super.addHandler(new ElementRecognizer.ByLocalName("seq"),
                GAMESeqHandler.GAME_SEQ_HANDLER_FACTORY);
        // <map_position>
//        super.addHandler(new ElementRecognizer.ByLocalName("map_position"),
//                GAMEMapPosHandler.GAME_MAP_POS_HANDLER_FACTORY);
        // <annotation>
        super.addHandler(new ElementRecognizer.ByLocalName("annotation"),
                GAMEAnnotationHandler.GAME_ANNOTATION_HANDLER_FACTORY);
        // <computational_analysis>
//        super.addHandler(new ElementRecognizer.ByLocalName("computational_analysis"),
//                GAMEComputationalHandler.GAME_COMP_HANDLER_FACTORY);
    }

    /**
     *  Description of the Method
     *
     *@param  nsURI             Description of the Parameter
     *@param  localName         Description of the Parameter
     *@param  qName             Description of the Parameter
     *@param  attrs             Description of the Parameter
     *@exception  SAXException  Description of the Exception
     */
    public void startElementHandler(
            String nsURI,
            String localName,
            String qName,
            Attributes attrs)
             throws SAXException 
    {
//        System.out.println("GAMEHandler startElementHandler called, localName: " + localName);
        // check the element type
        if (!(localName.equals("game"))) {
            throw new SAXException("first element of file is not a game element");
        }

        // check file version
        String version = attrs.getValue("version");
        if (!(version.equals("1.2"))) {
            throw new SAXException("GAME version is not 1.2");
        }

        // inform listener of new sequence
        try {
            listener.startSequence();
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEAnnotationHandler.");
        }
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler)
        throws SAXException
    {
        // indicate end of sequence
        try {
            listener.endSequence();
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("error in GAMEAnnotationHandler.");
        }
    }
}

