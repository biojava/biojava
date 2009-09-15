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
import java.util.Collection;
import java.util.List;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * seq_feature
 * 
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVESeqFeatureHandler
               extends StAXFeatureHandler implements AGAVEFeatureCallbackItf

{
  public static final StAXHandlerFactory AGAVE_SEQ_FEATURE_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVESeqFeatureHandler(staxenv);
    }
  };


 AGAVESeqFeatureHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("seq_feature", true);
    // setup handlers
       //
       super.addHandler(new ElementRecognizer.ByLocalName("classification"),
         AGAVEClassificationHandler.AGAVE_CLASSIFICATION_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("note"),
         AGAVENotePropHandler.AGAVE_NOTE_PROP_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("seq_location"),
         AGAVESeqLocationPropHandler.AGAVE_SEQ_LOCATION_PROP_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("xrefs"),
         AGAVEXrefsPropHandler.AGAVE_XREFS_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("evidence"),
         AGAVEEvidenceHandler.AGAVE_EVIDENCE_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("qualifier"),
         AGAVEQualifierPropHandler.AGAVE_QUALIFIER_PROP_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("seq_feature"),
         AGAVESeqFeatureHandler.AGAVE_SEQ_FEATURE_HANDLER_FACTORY);

      super.addHandler(new ElementRecognizer.ByLocalName("related_annot"),
         AGAVERelatedAnnotPropHandler.AGAVE_RELATED_ANNOT_PROP_HANDLER_FACTORY);

  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
                throws SAXException
  {
      try{
      featureListener.startFeature(featureTemplate);
      boolean forFeature = true ;
      setProperty( "element_id",  attrs.getValue("element_id"), forFeature ) ;
      setProperty( "feature_type",  attrs.getValue("feature_type"), forFeature ) ;
      setProperty( "label",  attrs.getValue("label"), forFeature ) ;

      featureTemplate.type = attrs.getValue("feature_type") ;
      }catch(Exception e)
      {
         throw new SAXException( e.getMessage() ) ;
      }
  }


  /*
   protected Feature.Template createTemplate() {
    // create Gene Template for this
    StrandedFeature.Template st = new StrandedFeature.Template();

    // assume feature set to describe a transcript
    st.type = "seq_feature";
    st.strand = StrandedFeature.UNKNOWN;
    // set up annotation bundle
    st.annotation = annot;
    st.location = new  Location.EmptyLocation();
    if( staxenv != null )
        staxenv. subFeatures .add( this ) ;

    return st;
  }*/


   public void addProperty(AGAVEProperty prop) {
      try{
	  String propType = prop.getPropType();
	  String value = prop.getValue();

	  // System.out.println(propType + " : " + value);

	  Object ob = null;
	  if (staxenv.featureTemplate.annotation.containsProperty(propType)) {
	      ob = staxenv.featureTemplate.annotation.getProperty(propType);
	  }

         if( ob != null && ob instanceof Collection) 
             ((Collection)ob).add( value ) ;
         else
         {
             List props = new ArrayList(1) ;
             props.add( value ) ;
             staxenv.featureTemplate.annotation.setProperty(propType, props)  ;
         }
      }catch (ChangeVetoException cve) {
          cve.printStackTrace() ;
      }
   }
  public void reportFeature(Location loc)
  {
      ((StrandedFeature.Template)featureTemplate).location = loc  ;
  }
  public void reportStrand(StrandedFeature.Strand strand)
  {
      ((StrandedFeature.Template)featureTemplate).strand = strand ;
  }

  public void addRelatedAnnot(AGAVERelatedAnnot prop)
  {
      try{
         Object ob =  UtilHelper.getProperty(staxenv.featureTemplate.annotation, "related_annot");
         if( ob != null )
             ((List)ob).add( prop ) ;
         else
         {
             List props = new ArrayList(1) ;
             props.add( prop ) ;
             staxenv.featureTemplate.annotation.setProperty("related_annot", props)  ;
         }
      }catch (ChangeVetoException cve) {
          cve.printStackTrace() ;
      }
   }


}
