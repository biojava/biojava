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
import java.util.ArrayList;
import java.util.List;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * transcript
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVETranscriptHandler
               extends StAXFeatureHandler implements AGAVEEvidenceCallbackItf
{

  public static final StAXHandlerFactory AGAVE_TRANSCRIPT_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVETranscriptHandler(staxenv);
    }
  };


 AGAVETranscriptHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("transcript", true);

    // setup handlers
       super.addHandler(new ElementRecognizer.ByLocalName("exons"),
         AGAVEExonsPropHandler.AGAVE_EXONS_PROP_HANDLER_FACTORY);
      super.addHandler(new ElementRecognizer.ByLocalName("cds"),
         AGAVECdsHandler.AGAVE_CDS_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("mrna"),
         AGAVEMrnaHandler.AGAVE_MRNA_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("predicted_protein"),
         AGAVEPredictedProteinHandler.AGAVE_PREDICTED_PROTEIN_HANDLER_FACTORY);
  }
  public void addElementId(String id)
  {
      try{
         Object ob = UtilHelper.getProperty(featureTemplate.annotation, "exons") ;
         if( ob == null ){
            ob = new ArrayList(1) ;
            featureTemplate.annotation.setProperty("exons", ob ) ;
         }
         ((List)ob).add( id) ;
      }catch(Exception e){
        e.printStackTrace() ;
      }
  }

   public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
      featureTemplate.type = "transcript" ;
  }


  /**
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "transcript";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = annot;
    st.location = new  Location.EmptyLocation();
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }**/


}
