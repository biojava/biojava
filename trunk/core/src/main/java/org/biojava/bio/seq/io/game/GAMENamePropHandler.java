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

import java.util.ListIterator;

import org.biojava.bio.seq.io.ParseException;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * StAX handler for the GAME &lt;name&gt; element.
 * derived from Thomas Down's PropDetailHandler
 *
 * @author David Huen
 * @author Thomas Down
 * @since 1.8
 */

public class GAMENamePropHandler extends StringElementHandlerBase {
    public static final StAXHandlerFactory GAME_NAME_PROP_HANDLER_FACTORY = new StAXHandlerFactory() {
	    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
//                System.out.println ("GAMENamePropHandler factory called.");
		return new GAMENamePropHandler(staxenv);
	    }
	} ;

  private StAXFeatureHandler staxenv;
  // this class is not derived from StAXPropertyHandler and doesn't inherit
  // any of the handlerStack maintenance code.  However it does offer
  // context either so it can be omitted from the stack.
  public GAMENamePropHandler(StAXFeatureHandler staxenv) {
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
//     System.out.println("GAMENamePropHandler.startElement entered.");
     super.startElement(nsURI, localName, qName, attrs, dm);
//     System.out.println("GAMENamePropHandler.startElement left.");
  }

  protected void setStringValue(String s)
        throws SAXException
  {
    int currLevel = staxenv.getLevel();
//    System.out.println("GAMENamePropHandler.setStringValue entered. currlevel: " + currLevel);

    if (currLevel >=1) {
      // search down stack for callback handler
      ListIterator li = staxenv.getHandlerStackIterator(currLevel);
//    System.out.println("GAMENamePropHandler.setStringValue entered. got ListIterator");      
      while (li.hasPrevious()) {
        Object ob = li.previous();
//    System.out.println("GAMENamePropHandler.setStringValue entered. got stack object");      
        if (ob instanceof GAMENameCallbackItf) {
          // we have a nesting handler, use it
//    System.out.println("GAMENamePropHandler.setStringValue calling back");              
          ((GAMENameCallbackItf) ob).NameSetStringValue(s);
          return;
        }
      }
    }
    //  default is to just set the name property to the string.
    try {
      staxenv.getFeatureListener().addFeatureProperty("name", s);
    }
    catch (ParseException pe) {
      throw new SAXException(pe);
    }
  }
}
