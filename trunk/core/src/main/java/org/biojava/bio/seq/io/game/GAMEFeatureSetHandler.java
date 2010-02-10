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
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * Handles the &lt;feature_set&gt; element
 *
 * @author David Huen
 * @since 1.2
 */
public class GAMEFeatureSetHandler
               extends StAXFeatureHandler
               implements GAMENameCallbackItf, GAMETranscriptCallbackItf  
{
  // <feature_span> is one of the worst elements in GAME.  The type of element
  // is only known from reading the nested <type> element.  The bloody
  // strand and coordinates are nested 4-deep within the element.
  // this element is implemented as a transcript.

  private Location transcriptRange;
  private StAXFeatureHandler staxenv;

  public static final StAXHandlerFactory GAME_FEATURESET_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEFeatureSetHandler(staxenv);
    }
  };

  GAMEFeatureSetHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    this.staxenv = staxenv;
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("feature_set", true);

    // setup handlers
       // <name>
       super.addHandler(new ElementRecognizer.ByLocalName("name"),
         GAMENamePropHandler.GAME_NAME_PROP_HANDLER_FACTORY);
       // <feature_span>
       super.addHandler(new ElementRecognizer.ByLocalName("feature_span"),
         GAMEFeatureSpanHandler.GAME_FEATURESPAN_HANDLER_FACTORY);
  }


  public void NameSetStringValue(String s) {
    // the id overrides name element by default.
    if (!featureTemplate.annotation.containsProperty("id")) {
      try {
       featureTemplate.annotation.setProperty("id", s);
      }
      catch (ChangeVetoException cve) {
        // baulk and discard exception
        System.out.println("GAMEFeatureSetHandler: change vetoed");
      }
    }
  }

  public void reportExon(RangeLocation range, StrandedFeature.Strand strand) {
    // set strand of transcript
    ((StrandedFeature.Template) featureTemplate).strand = strand;

    // check extend current range with new exon
    if (transcriptRange == null)
      transcriptRange = range;
    else
      transcriptRange = transcriptRange.union(range);
  }

  protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "transcript";
    st.strand = StrandedFeature.UNKNOWN;

    // set up annotation bundle
    st.annotation = new SmallAnnotation();

    return st;
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
	 throws SAXException
  {
    // handle the attributes
    String s = (String) attrs.getValue("id");
    if (s != null) {
      try {
        featureTemplate.annotation.setProperty("id", s);
      }
      catch (ChangeVetoException cve) {
        // just ignore change
        System.err.println("GAMEFeatureSetHandler: change vetoed");
      }
    }
  }

  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
  {
    // set the extent of the transcript to its limits
    if (transcriptRange != null) {
//      featureTemplate.location = new RangeLocation(transcriptRange.getMin(), transcriptRange.getMax());
      featureTemplate.location = transcriptRange;      
    }
    else
      featureTemplate.location = Location.empty;

    // tell nesting feature about my extent and strand.
    int currLevel = staxenv.getLevel();
//    System.out.println("GAMEFeatureSetHandler.endElement entered. currlevel: " + currLevel);
 
    if (currLevel >=1) {
      // search down stack for callback handler
      ListIterator li = staxenv.getHandlerStackIterator(currLevel);
//    System.out.println("GAMEFeatureSetHandler.endElement entered. got ListIterator");
      while (li.hasPrevious()) {
        Object ob = li.previous();
//    System.out.println("GAMEFeatureSetHandler.endElement entered. got stack object");
        if (ob instanceof GAMEFeatureCallbackItf) {
          // we have a nesting handler, use it
//    System.out.println("GAMEFeatureSetHandler.endElement calling back");
          ((GAMEFeatureCallbackItf) ob).reportFeature(featureTemplate.location);
          ((GAMEFeatureCallbackItf) ob).reportStrand(((StrandedFeature.Template) featureTemplate).strand);
          return;
        }
      }
    }

  }
}

