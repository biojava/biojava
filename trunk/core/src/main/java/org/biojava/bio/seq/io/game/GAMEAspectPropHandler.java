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
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handles the GAME &lt;aspect&gt; element
 *
 * @author David Huen
 * @since 1.8
 */
public class GAMEAspectPropHandler extends StAXPropertyHandler {
  // the <seq> element supplies clone name and length.
  // other data includes a description of the sequence.
  // we will stuff the name as clone_name in an annotation.

  // set up factory method
  public static final StAXHandlerFactory GAME_ASPECT_PROP_HANDLER_FACTORY 
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEAspectPropHandler(staxenv);
    }
  };

  GAMEAspectPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);

    // setup handlers
    super.addHandler(new ElementRecognizer.ByLocalName("name"),
      GAMENamePropHandler.GAME_NAME_PROP_HANDLER_FACTORY);
  }

  public void startElement(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs,
                DelegationManager dm)
	 throws SAXException
  {
    // all seems in order, go to superclass method
    super.startElement(nsURI, localName, qName, attrs, dm);
  }
}

