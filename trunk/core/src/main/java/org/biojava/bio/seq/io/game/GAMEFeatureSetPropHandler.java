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

/**
 * Handles the GAME <feature_set> element
 *
 * @author David Huen
 * @since 1.8
 */
public class GAMEFeatureSetPropHandler extends StAXPropertyHandler {
  // while <feature_set> doesn't invoke creation of features, elements nested
  // by it do.  This element just acts as a wrapper.
  // there are attributes on this element.

  // set up factory method
  public static final StAXHandlerFactory GAME_FEATURESET_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEFeatureSetPropHandler(staxenv);
    }
  };

  GAMEFeatureSetPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("feature_set", false);

    // setup handlers
    // rather daft one, this: there's already an id attribute.
    super.addHandler(new ElementRecognizer.ByLocalName("name"),
      GAMENamePropHandler.GAME_NAME_PROP_HANDLER_FACTORY);
    // rather daft one, this: there's already an id attribute.
    super.addHandler(new ElementRecognizer.ByLocalName("feature_span"),
      GAMEFeatureSpanHandler.GAME_FEATURESPAN_HANDLER_FACTORY);
  }

/*
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
	 throws SAXException
  {
  }
*/
}

