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

import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;

/**
 * Handles the GAME &lt;annotation&gt; element
 *
 * @author David Huen
 * @since 1.2
 */
public class GAMEAnnotationHandler 
               extends StAXFeatureHandler 
               implements GAMEFeatureCallbackItf {
  // all Gadfly data concerning a single "gene" appears in
  // single <annotation>.
  // set up factory method
  public static final StAXHandlerFactory GAME_ANNOTATION_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new GAMEAnnotationHandler(staxenv);
    }
  };

  private Location range = Location.empty;

  GAMEAnnotationHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("annotation", true);

    // setup handlers
       // <seq>: never seen it used yet.
//       super.addHandler(new ElementRecognizer.ByLocalName("seq"),
//         GAMESeqPropHandler.GAME_SEQ_PROP_HANDLER_FACTORY);
       // <map_position>
//       super.addHandler(new ElementRecognizer.ByLocalName("map_position"),
//         GAMEMapPosPropHandler.GAME_MAP_POS_PROP_HANDLER_FACTORY);
       // <gene>
       super.addHandler(new ElementRecognizer.ByLocalName("gene"),
         GAMEGenePropHandler.GAME_GENE_PROP_HANDLER_FACTORY);
       // <feature_set>
       super.addHandler(new ElementRecognizer.ByLocalName("feature_set"),
         GAMEFeatureSetHandler.GAME_FEATURESET_HANDLER_FACTORY);
       // <dbxref>
       super.addHandler(new ElementRecognizer.ByLocalName("dbxref"),
         GAMEDbxrefPropHandler.GAME_DBXREF_PROP_HANDLER_FACTORY);
       // <Aspect>
//       super.addHandler(new ElementRecognizer.ByLocalName("aspect"),
//         GAMEAspectPropHandler.GAME_ASPECT_PROP_HANDLER_FACTORY);
  }

  protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template gt = new StrandedFeature.Template();

    // set up annotation bundle
    gt.type = "gene";
    gt.source = "";
    gt.location = Location.empty;
    gt.annotation = new SmallAnnotation();
    gt.strand = StrandedFeature.UNKNOWN;

    return gt;
  }

  public void reportFeature(Location loc)
  {
//    System.out.println("GAMEAnnotationHandler location is " + loc);
    // accumulate locations of features here.
    range = range.union(loc);
//    System.out.println("GAMEAnnotationHandler after union is  " + range);
  }

  public void reportStrand(StrandedFeature.Strand strand)
  {
    // obtains strand from elements that are in the know.
    ((StrandedFeature.Template) featureTemplate).strand = strand;
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
  {
    String annotationId =  attrs.getValue("id");
    if (annotationId != null) {
      // stuff Gadfly annotation id into our annotation bundle for info.
//      System.out.println("GAMEAnnotationHandler is setting id to " + annotationId);
      try {
         featureTemplate.annotation.setProperty(
                          "annotation_id", annotationId);
      }
      catch (ChangeVetoException cae) {
        System.err.println("GAMEAnnotationHandler: veto exception caught.");
      }
    }
//    System.out.println("GAMEAnnotationHandler.startElementHandler: leaving.");
  }

  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
  {
    // finalise the sequence extent to encompass all reported features
    if (range != Location.empty)
      featureTemplate.location = new RangeLocation(
                                       range.getMin(), 
                                       range.getMax());
  }

}

