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

import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handles the GAME &lt;map_position&gt; element
 * Currently, it just ignores it!
 *
 * @author David Huen
 * @since 1.8
 */
public class GAMEMapPosPropHandler extends StAXPropertyHandler {
  // the <map_position> element supplies details of the map
  // position of the clone.  It can be the position on
  // some canonical master sequence(tile) or cytology(cytological) depending
  // on the type attribute.

  // set up factory method
  public static final StAXHandlerFactory GAME_MAP_POS_PROP_HANDLER_FACTORY 
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEMapPosPropHandler(staxenv);
    }
  };

  GAMEMapPosPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);

    setHandlerCharacteristics("map_position", false);
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
	 throws SAXException
  {
  }
}

