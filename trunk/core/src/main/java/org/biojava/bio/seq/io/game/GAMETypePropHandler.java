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

package org.biojava.bio.seq.io.game;

import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for GAME &lt;type&gt; elements.
 * derived from Thomas Down's PropDetailHandler
 *
 * @author David Huen
 * @author Thomas Down
 * @since 1.8
 */
public class GAMETypePropHandler extends StringElementHandlerBase {
  // The <type> handler is context sensitive as the meaning of type depends on
  // the element in which it is nested.
  public static final StAXHandlerFactory GAME_TYPE_PROP_HANDLER_FACTORY = new StAXHandlerFactory() {
	  public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                 return new GAMETypePropHandler(staxenv);
	  }
            } ;

  private StAXFeatureHandler staxenv;

  public GAMETypePropHandler(StAXFeatureHandler staxenv) {
    super();
    this.staxenv = staxenv;
  }

  public void startElement(String nsURI,
                                    String localName,
                                    String qName,
                                    Attributes attrs,
                                    DelegationManager dm)
                     throws SAXException
  {

    super.startElement(nsURI, localName, qName, attrs, dm);
  }

  protected void setStringValue(String s)
        throws SAXException
  {
      // if this is <feature_span, the string will decide the type of feature to create
      String trimmed = s.trim();

      if (trimmed.equals("translate offset")) {
        // create a start codon annotation
//        System.out.println("setting ATG");
        staxenv.featureTemplate.type = "ATG";
      }
      else if (trimmed.equals("exon")) {
        //  feature is an exon
//        System.out.println("Setting exon");
        staxenv.featureTemplate.type = "exon";
      }
   }
}
