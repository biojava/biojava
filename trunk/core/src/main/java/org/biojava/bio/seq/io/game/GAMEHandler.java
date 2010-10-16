/**
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

package org.biojava.bio.seq.io.game;

import org.biojava.utils.stax.DelegationManager;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handles the root GAME element
 *
 * @author David Huen
 * @since 1.8
 */
public class GAMEHandler extends StAXFeatureHandler {
  // there is only one GAME element encompassing the entire file
  // in Gadfly GAME files.  We create a single feature template
  // collect all its info into one annotation bundle.

  public GAMEHandler() {
    // set up environment
    setHandlerCharacteristics("game", true);

    // setup handlers
       // <seq>
       super.addHandler(new ElementRecognizer.ByLocalName("seq"), 
         GAMESeqPropHandler.GAME_SEQ_PROP_HANDLER_FACTORY); 
       // <map_position>
       super.addHandler(new ElementRecognizer.ByLocalName("map_position"), 
         GAMEMapPosPropHandler.GAME_MAP_POS_PROP_HANDLER_FACTORY); 
       // <annotation>
       super.addHandler(new ElementRecognizer.ByLocalName("annotation"), 
         GAMEAnnotationHandler.GAME_ANNOTATION_HANDLER_FACTORY); 
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs,
                DelegationManager dm)
	 throws SAXException
  {
    // check the element type
    if (!(localName.equals( "game")) )
      throw new SAXException("first element of file is not a game element");

    // check file version
    String version = attrs.getValue("version");
    if (!(version.equals("1.001")) )
      throw new SAXException("GAME version is not 1.001");

//    System.out.println("GAMEHandler startElementHandler called, localName: " + localName);
  }
}

