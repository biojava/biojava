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
package org.biojava.bio.seq.io.agave;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

 /**
  * sequence_map
  *
  * @author Hanning Ni    Doubletwist Inc
  */
public class AGAVESeqMapHandler
               extends StAXFeatureHandler

{
  public static final StAXHandlerFactory AGAVE_SEQ_MAP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVESeqMapHandler(staxenv);
    }
  };


  AGAVESeqMapHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("sequence_map", true);

    // setup handlers
       //
       super.addHandler(new ElementRecognizer.ByLocalName("note"),
         AGAVENotePropHandler.AGAVE_NOTE_PROP_HANDLER_FACTORY);
       //
     //  super.addHandler(new ElementRecognizer.ByLocalName("computation"),
     //    AGAVEComputationHandler.AGAVE_COMPUTATION_HANDLER_FACTORY);
      super.addHandler(new ElementRecognizer.ByLocalName("annotations"),
         AGAVEAnnotationsHandler.AGAVE_ANNOTATIONS_HANDLER_FACTORY);

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
      setProperty( "label",  attrs.getValue("label") , forFeature) ;
      }catch(Exception e){
         throw new SAXException( e.getMessage() ) ;
      }
  }

  /**
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "sequence_map";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = annot;
    st.location = new  Location.EmptyLocation();
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }**/



}

