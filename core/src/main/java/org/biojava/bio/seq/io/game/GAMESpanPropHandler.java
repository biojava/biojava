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

import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.PointLocation;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.stax.IntElementHandlerBase;
import org.biojava.utils.stax.StAXContentHandler;

/**
 * Handles the GAME &lt;span&gt; element
 * Currently, it just ignores it!
 *
 * @author David Huen
 * @since 1.2
 */
public class GAMESpanPropHandler 
               extends StAXPropertyHandler {
  // the <span> element supplies limits of a sequence span.
  // unfortunately, the spans can be either numeric or
  // alphanumeric (with cytological map_position).
  // set up factory method
  public static final StAXHandlerFactory GAME_SPAN_PROP_HANDLER_FACTORY 
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMESpanPropHandler(staxenv);
    }
  };

  private int start = 0;
  private int stop = 0;
  private StAXFeatureHandler staxenv;

  GAMESpanPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("span", true);
   
    // cache environment: this is of PREVIOUS feature handler as
    // delegation is invoked in StaxFeatureHandler itself which means
    // that the this pointer that is passed is the Feature one.
    this.staxenv = staxenv;

    // setup handlers
    super.addHandler(new ElementRecognizer.ByLocalName("start"),
//      GAMEStartEndPropHandler.GAME_STARTEND_PROP_HANDLER_FACTORY);
      new StAXHandlerFactory() {
           public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
             return new StartHandler(); }
      }
    );

    super.addHandler(new ElementRecognizer.ByLocalName("end"),
//      GAMEStartEndPropHandler.GAME_STARTEND_PROP_HANDLER_FACTORY);
      new StAXHandlerFactory() {
           public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
             return new StopHandler(); }
      }
    );
  }

  private class StartHandler extends IntElementHandlerBase
  {
    protected void setIntValue(int val)
    {
      start = val;
    }
  }
 
  private class StopHandler extends IntElementHandlerBase
  {
    protected void setIntValue(int val)
    {
      stop = val;
    }
  }

/*
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
	 throws SAXException
  {
    System.out.println("GAMESpanPropHandler.startElementHandler entered.");
  }
*/
  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
  {
    // check that it IS a StrandedFeature.Template
    boolean isStrandedTemplate = staxenv.featureTemplate instanceof StrandedFeature.Template;

    Feature.Template templ = staxenv.featureTemplate;

    // go set strandedness and range
    if (start < stop) {
      templ.location = new RangeLocation(start + 1, stop);
      if (isStrandedTemplate)
        ((StrandedFeature.Template) templ).strand = StrandedFeature.POSITIVE;
    }
    else if (start > stop) {
      staxenv.featureTemplate.location = new RangeLocation(stop + 1, start);
      if (isStrandedTemplate)
        ((StrandedFeature.Template) templ).strand = StrandedFeature.NEGATIVE;
    }
    else {
      staxenv.featureTemplate.location = new PointLocation(start);
      if (isStrandedTemplate)
        ((StrandedFeature.Template) templ).strand = StrandedFeature.UNKNOWN;
    } 
  }
}

