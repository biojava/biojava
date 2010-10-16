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
import java.util.ListIterator;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 *
 *
 */
public class AGAVERelatedAnnotPropHandler extends StAXPropertyHandler implements AGAVEEvidenceCallbackItf{


   public static final StAXHandlerFactory AGAVE_RELATED_ANNOT_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVERelatedAnnotPropHandler(staxenv);
    }
   };
   private AGAVERelatedAnnot related_annot ;
   AGAVERelatedAnnotPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("related_annot", true);

    super.addHandler(new ElementRecognizer.ByLocalName("element_id"),
         AGAVEElementIdPropHandler.AGAVE_ELEMENT_ID_PROP_HANDLER_FACTORY);
    super.addHandler(new ElementRecognizer.ByLocalName("sci_property"),
         AGAVESciPropertyPropHandler.AGAVE_SCI_PROPERTY_PROP_HANDLER_FACTORY);
    related_annot = new AGAVERelatedAnnot() ;
  }
  public void addElementId(String id)
  {
      related_annot.addElementId(id) ;
  }
  public void addProperty(AGAVEProperty prop)
  {
      related_annot.addProp( prop ) ;
  }
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
      related_annot.setScore(attrs.getValue( "score" )  );
      related_annot.setRel( attrs.getValue( "rel" )   )  ;
  }


  public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
       int currLevel = staxenv.getLevel();
       if (currLevel >=1) {
           ListIterator li = staxenv.getHandlerStackIterator(currLevel);
           while( li.hasPrevious() )
          {
              Object ob =   li.previous() ;
              if(  ( ob instanceof AGAVECompResultHandler ) ||
                  ( ob instanceof AGAVESeqFeatureHandler ) )
              {
                  ( (AGAVECompResultHandler) ob ).addRelatedAnnot( related_annot ) ;
                   return ;
              }
           }
       }

  }


}
