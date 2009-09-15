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
 *
 * @author Hanning Ni     Doubletwist Inc
 */
public class AGAVEClassificationHandler
               extends StAXFeatureHandler
               implements AGAVEIdAliasCallbackItf
{
  public static final StAXHandlerFactory AGAVE_CLASSIFICATION_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEClassificationHandler(staxenv);
    }
  };

 private List idAliases ;
 AGAVEClassificationHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("classification", true);

    // setup handlers
       //
       super.addHandler(new ElementRecognizer.ByLocalName("description"),
         AGAVEDescPropHandler.AGAVE_DESC_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("id_alias"),
         AGAVEIdAliasPropHandler.AGAVE_ID_ALIAS_PROP_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("evidence"),
         AGAVEEvidenceHandler.AGAVE_EVIDENCE_HANDLER_FACTORY);


  }

   public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
                throws SAXException
  {
      try{
      featureListener.startFeature( featureTemplate );
      boolean forFeature = true ;
      setProperty( "system",  attrs.getValue("system"), forFeature ) ;
      setProperty( "id",  attrs.getValue("id"), forFeature ) ;
      setProperty( "type",  attrs.getValue("type") , forFeature) ;
      setProperty( "assigned_by",  attrs.getValue("assigned_by"), forFeature ) ;
      }catch(Exception e){
         throw new SAXException( e.getMessage() ) ;
      }
  }


  public void addIdAlias(AGAVEIdAlias id)
  {
      if( idAliases == null )
          idAliases = new ArrayList(1) ;
      idAliases.add( id ) ;
  }
  /*
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "classification";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = annot;
    st.location = new  Location.EmptyLocation();
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }*/

  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
              throws SAXException
  {
      try{
          if( idAliases != null )
              annot.setProperty("id_alias", idAliases) ;
      }catch(Exception e){
        throw new SAXException( e.getMessage() ) ;
      }
  }


}
