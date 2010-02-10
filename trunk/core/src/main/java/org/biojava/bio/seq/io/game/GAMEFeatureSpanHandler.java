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

import java.util.ListIterator;

import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;

/**
 * Handles the &lt;feature_span&gt; element
 *
 * @author David Huen
 * @since 1.8
 */
public class GAMEFeatureSpanHandler extends StAXFeatureHandler {
  // <feature_span> is one of the worst elements in GAME.  The type of element
  // is only known from reading the nested <type> element.  The bloody
  // strand and coordinates are nested 4-deep within the element.
  // this element could be a start codon feature or an exon.
  // Also, the lower value coordinate is one less than than the equivalent
  // EMBL coordinate but the higher value is the same.
  private StAXFeatureHandler staxenv;

  public static final StAXHandlerFactory GAME_FEATURESPAN_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEFeatureSpanHandler(staxenv);
    }
  };

  GAMEFeatureSpanHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    this.staxenv = staxenv;
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("feature_span", false);

    // setup handlers
       // <type>
       super.addHandler(new ElementRecognizer.ByLocalName("type"),
         GAMETypePropHandler.GAME_TYPE_PROP_HANDLER_FACTORY);
       // <seq_relationship>
       super.addHandler(new ElementRecognizer.ByLocalName("seq_relationship"),
         GAMESeqRelPropHandler.GAME_SEQREL_PROP_HANDLER_FACTORY);
  }

  protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template ft = new StrandedFeature.Template();

    // set up annotation bundle
    ft.annotation = new SmallAnnotation();
    ft.strand = StrandedFeature.UNKNOWN;

    return ft;
  }


  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
  {
//    System.out.println("GAMEFeatureSpanHandler.startElementHandler entered.");

      // pick up id and save it in annotation bundle
      String featureId = attrs.getValue("id");

      try {
        if (featureId != null)
          featureTemplate.annotation.setProperty("id", featureId);
      }
      catch (ChangeVetoException cve) {
        System.err.println("GAMEFeatureSpanHandler. Change blocked.");
      }
  }

  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
  {
    // we only do this for exons
    if (!( ((String) featureTemplate.type).equals("exon") )) return;

    // update transcript limits
    // get iterator to callbackStack of PREVIOUS FeatureHandler
    int currLevel = staxenv.getLevel();
 
    if (currLevel >=1) {
      // search down stack for callback handler
      ListIterator li = staxenv.getHandlerStackIterator(currLevel);

      while (li.hasPrevious()) {
        Object ob = li.previous();
        if (ob instanceof GAMETranscriptCallbackItf) {
          // we have a nesting handler, use it
          ((GAMETranscriptCallbackItf) ob).reportExon(
              (RangeLocation) ((StrandedFeature.Template) featureTemplate).location,
              ((StrandedFeature.Template) featureTemplate).strand);
          return;
        }
      }
    }
    
  }

}

