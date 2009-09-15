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
public class AGAVEFragmentOrientationHandler
               extends StAXFeatureHandler implements AGAVEFeatureCallbackItf

{
  public static final StAXHandlerFactory AGAVE_FRAGMENT_ORIENTATION_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEFragmentOrientationHandler(staxenv);
    }
  };


  AGAVEFragmentOrientationHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("fragment_orientation", true);

    // setup handlers
        //
       super.addHandler(new ElementRecognizer.ByLocalName("bio_sequence"),
         AGAVEBioSeqHandler.AGAVE_BIO_SEQ_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("map_location"),
         AGAVEMapLocationPropHandler.AGAVE_MAP_LOCATION_PROP_HANDLER_FACTORY);

  }


  public void reportFeature(Location loc)
  {
    ((StrandedFeature.Template) featureTemplate).location = loc;
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
                throws SAXException
  {
     try{
      featureListener.startFeature( featureTemplate);
      String strand = attrs.getValue("on_complement_strand");
      boolean forFeature = true ;
      setProperty( "on_complement_strand",  strand , forFeature) ;
      setProperty( "has_5prime_end",  attrs.getValue("has_5prime_end") , forFeature) ;
      setProperty( "has_3prime_end",  attrs.getValue("has_3prime_end"), forFeature ) ;
      setProperty( "is_all_BAC_vect",  attrs.getValue("is_all_BAC_vect") , forFeature) ;
      if( strand.equalsIgnoreCase("true") )
        ((StrandedFeature.Template) featureTemplate).strand =  StrandedFeature.NEGATIVE ;
      else if( strand.equalsIgnoreCase("false") )
         ((StrandedFeature.Template) featureTemplate).strand =  StrandedFeature.POSITIVE ;
      else
        ((StrandedFeature.Template) featureTemplate).strand =  StrandedFeature.UNKNOWN ;
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
