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
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *
 * @author Hanning Ni     Doubletwist Inc
 */
public class AGAVECompResultHandler
               extends StAXFeatureHandler implements AGAVEFeatureCallbackItf

{
  public static final StAXHandlerFactory AGAVE_COMP_RESULT_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVECompResultHandler(staxenv);
    }
  };


  AGAVECompResultHandler(StAXFeatureHandler staxenv) {
    // setup up environment stuff
    super( staxenv );
    featureListener = staxenv.featureListener;
    setHandlerCharacteristics("comp_result", true);

    // setup handlers
        //
       super.addHandler(new ElementRecognizer.ByLocalName("note"),
         AGAVENotePropHandler.AGAVE_NOTE_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("match_desc"),
         AGAVEMatchDescPropHandler.AGAVE_MATCH_DESC_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("match_align"),
         AGAVEMatchAlignPropHandler.AGAVE_MATCH_ALIGN_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("query_region"),
         AGAVEQueryRegionPropHandler.AGAVE_QUERY_REGION_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("match_region"),
         AGAVEMatchRegionPropHandler.AGAVE_MATCH_REGION_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("result_property"),
         AGAVEResultPropertyPropHandler.AGAVE_RESULT_PROPERTY_PROP_HANDLER_FACTORY);
       //
       super.addHandler(new ElementRecognizer.ByLocalName("result_group"),
         AGAVEResultGroupHandler.AGAVE_RESULT_GROUP_HANDLER_FACTORY);
       //
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
      featureListener.startFeature( featureTemplate );
      boolean forFeature = true ;
      setProperty( "element_id",  attrs.getValue("element_id") , forFeature) ;
      setProperty( "result_id",  attrs.getValue("result_id") , forFeature) ;
      setProperty( "group_order",  attrs.getValue("group_order") , forFeature) ;
      setProperty("result_type",  attrs.getValue("result_type"), forFeature ) ;
      setProperty( "feature_type",  attrs.getValue("feature_type"), forFeature ) ;
      setProperty( "on_complement_strand",  attrs.getValue("on_complement_strand") , forFeature) ;
      setProperty( "confidence",  attrs.getValue("confidence"), forFeature ) ;
      setProperty( "align_length",  attrs.getValue("align_length") , forFeature) ;
      setProperty( "align_unit",  attrs.getValue("align_unit") , forFeature) ;
      String strand = attrs.getValue("on_complement_strand") ;
      if( strand.equalsIgnoreCase("true") )
        ((StrandedFeature.Template) featureTemplate).strand =  StrandedFeature.NEGATIVE ;
      else if( strand.equalsIgnoreCase("false") )
         ((StrandedFeature.Template) featureTemplate).strand =  StrandedFeature.POSITIVE ;
      else
        ((StrandedFeature.Template) featureTemplate).strand =  StrandedFeature.UNKNOWN ;
      featureTemplate.type = attrs.getValue("result_type") ;
    }catch(Exception e){
        throw new SAXException( e.getMessage() ) ;
    }
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
  public void addProperty(AGAVEProperty prop)
  {
      try{
         Object ob = UtilHelper.getProperty(staxenv.featureTemplate.annotation,"result_property");
         if( ob != null )
             ((List)ob).add( prop ) ;
         else
         {
             List props = new ArrayList(1) ;
             props.add( prop ) ;
             staxenv.featureTemplate.annotation.setProperty("result_property", props)  ;
         }
      }catch (ChangeVetoException cve) {
         cve.printStackTrace() ;
      }
   }

  public void addRelatedAnnot(AGAVERelatedAnnot prop)
  {
      try{
         Object ob = UtilHelper.getProperty(staxenv.featureTemplate.annotation, "related_annot");
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

