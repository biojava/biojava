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

import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handles the GAME &lt;seq&gt; element
 *
 * @since David Huen
 * @since 1.2
 */
public class GAMESeqPropHandler 
               extends StAXPropertyHandler 
               implements GAMENameCallbackItf {
  // the <seq> element supplies clone name and length.
  // other data includes a description of the sequence.
  // we will stuff the name as clone_name in an annotation.

  // set up factory method
  public static final StAXHandlerFactory GAME_SEQ_PROP_HANDLER_FACTORY 
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
//      System.out.println("GAMESeqPropHandler factory called.");
      StAXContentHandler temp = new GAMESeqPropHandler(staxenv);
//      System.out.println("GAMESeqPropHandler factory created " + temp);  
//      System.out.println("");
//      if (temp == null) System.out.println("GAMESeqPropHandler instantiation failed");
//      return new GAMESeqPropHandler(staxenv);
      return temp;
    }
  };

  GAMESeqPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
//    System.out.println("GAMESeqPropHandler constructor called.");
    setHandlerCharacteristics("seq", true);

    // setup handlers
    // <name>
    super.addHandler(new ElementRecognizer.ByLocalName("name"),
      GAMENamePropHandler.GAME_NAME_PROP_HANDLER_FACTORY);
    // <description>
    super.addHandler(new ElementRecognizer.ByLocalName("description"),
      GAMEDescriptionPropHandler.GAME_DESCRIPTION_PROP_HANDLER_FACTORY);
    // <residues>
    super.addHandler(new ElementRecognizer.ByLocalName("residues"),
      GAMEResiduesPropHandler.GAME_RESIDUES_PROP_HANDLER_FACTORY);
//    System.out.println("GAMESeqPropHandler constructor: leaving now.");
  }

  public void NameSetStringValue(String s) {
    if (!staxenv.featureTemplate.annotation.containsProperty("id")) {
      // set gene name
      try {
       staxenv.featureTemplate.annotation.setProperty("id", s.trim());
      }
      catch (ChangeVetoException cve) {
        // baulk and discard exception
        System.err.println("GAMEGenPropHandler: change vetoed");
      }
    }
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
	 throws SAXException
  {
    // pick up the sequence details set them on the sequence object
    // I will assume the length is equivalent to the coordinate range
    String lengthString =  attrs.getValue("length");

    if (lengthString != null) {
      staxenv.featureTemplate.location = new RangeLocation(1, Integer.parseInt(lengthString));
    }
  }
}

