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

import org.xml.sax.SAXException;

 /**
  *
  * @author Hanning Ni    Doubletwist Inc
  */
public class AGAVEEvidenceHandler
               extends StAXFeatureHandler implements AGAVEEvidenceCallbackItf

{
  public static final StAXHandlerFactory AGAVE_EVIDENCE_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEEvidenceHandler(staxenv);
    }
  };
 private List element_ids ;

 AGAVEEvidenceHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("evidence", true);

    // setup handlers
       //
       super.addHandler(new ElementRecognizer.ByLocalName("element_id"),
         AGAVEElementIdPropHandler.AGAVE_ELEMENT_ID_PROP_HANDLER_FACTORY);
       super.addHandler(new ElementRecognizer.ByLocalName("comp_result"),
         AGAVECompResultHandler.AGAVE_COMP_RESULT_HANDLER_FACTORY);

  }

  public void addElementId(String id)
  {
      if( element_ids == null )
         element_ids = new ArrayList(1) ;
      element_ids.add( id ) ;
  }
  /*
  protected Feature.Template createTemplate() {
    Feature.Template st = new Feature.Template();
    st.type = "evidence";
    st.annotation = annot;
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  } */
  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
              throws SAXException
  {
      try{
          if( element_ids != null )
              annot.setProperty("element_ids", element_ids ) ;
      }catch(Exception e){
          throw new SAXException( e.getMessage() ) ;
      }
  }
}

