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
package org.biojava.bio.seq.io.agave;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEFragmentOrderHandler
               extends StAXFeatureHandler

{
  public static final StAXHandlerFactory AGAVE_FRAGMENT_ORDER_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEFragmentOrderHandler(staxenv);
    }
  };


  AGAVEFragmentOrderHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("fragment_order", true);

    // setup handlers
       //
       super.addHandler(new ElementRecognizer.ByLocalName("fragment_orientation"),
         AGAVEFragmentOrientationHandler.AGAVE_FRAGMENT_ORIENTATION_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("bio_sequence"),
         AGAVEBioSequenceHandler.AGAVE_BIO_SEQUENCE_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("map_location"),
         AGAVEMapLocationPropHandler.AGAVE_MAP_LOCATION_PROP_HANDLER_FACTORY);

  }

  public void reportStrand(StrandedFeature.Strand strand)
  {
    // obtains strand from elements that are in the know.
    ((StrandedFeature.Template) featureTemplate).strand = strand;
  }
  public void reportFeature(Location loc)
  {
    ((StrandedFeature.Template) featureTemplate).location = loc;
  }



  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
                throws SAXException
  {
      try{
      featureListener.startSequence();
      boolean forFeature = false ;
      setProperty( "group_id",  attrs.getValue("group_id") , forFeature) ;
      setProperty( "length",  attrs.getValue("length") , forFeature) ;
      }catch(Exception e){
          throw new SAXException( e.getMessage() ) ;
      }
  }

  /**
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "fragment_order";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = annot;
    st.location = new  Location.EmptyLocation();
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }**/



}

