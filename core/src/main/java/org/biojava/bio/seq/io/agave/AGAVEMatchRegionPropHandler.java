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
 * Deals with match_region
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEMatchRegionPropHandler
               extends StAXPropertyHandler implements AGAVEDbIdCallbackItf, AGAVEEvidenceCallbackItf
{
  // set up factory method
  public static final StAXHandlerFactory AGAVE_MATCH_REGION_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEMatchRegionPropHandler(staxenv);
    }
  };
  private AGAVEMatchRegion match_region = new AGAVEMatchRegion() ;

  AGAVEMatchRegionPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("match_region", true);

    super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);
    super.addHandler(new ElementRecognizer.ByLocalName("element_id"),
         AGAVEElementIdPropHandler.AGAVE_ELEMENT_ID_PROP_HANDLER_FACTORY);
    //bio_sequence
    //super.addHandler(new ElementRecognizer.ByLocalName("bio_sequence"),
      //   AGAVEBioSequencePropHandler.AGAVE_BIO_SEQUENCE_PROP_HANDLER_FACTORY);


  }
  public void addDbId( AGAVEDbId id)
  {
      match_region.setDbId(id) ;
  }
  public void addElementId(String id)
  {
      match_region.setElementId( id ) ;
  }
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
     match_region.setStart( new Integer( attrs.getValue("start") ) .intValue()  );
     match_region.setEnd( new Integer( attrs.getValue("end") ) .intValue()  );
  }

   public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
      try{
          staxenv.featureTemplate.annotation.setProperty("match_region", match_region) ;
      }catch(Exception e){
          throw new SAXException( e.getMessage() ) ;
      }
  }

  }
