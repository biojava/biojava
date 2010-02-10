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

import org.biojava.bio.symbol.RangeLocation;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * 
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEQueryRegionPropHandler extends StAXFeatureHandler
                                         implements AGAVEDbIdCallbackItf {

   // set up factory metho
  public static final StAXHandlerFactory AGAVE_QUERY_REGION_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEQueryRegionPropHandler(staxenv);
    }
  };
   private AGAVEQueryRegion region ;
   AGAVEQueryRegionPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("query_region", true);
    super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);
    region = new AGAVEQueryRegion() ;
   }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
      region.setStart( new Integer(attrs.getValue( "start" )).intValue() );
      region.setEnd( new Integer(attrs.getValue( "end" ) ).intValue() );
  }

  public void addDbId(AGAVEDbId db_id)
  {
     try{
       region.setDbId( db_id ) ;
     }catch(Exception e){
       e.printStackTrace() ;
     }
  }


   public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
        try{
           staxenv.featureTemplate.annotation.setProperty("query_region", region) ;
        }catch(Exception e){
          throw new SAXException( e.getMessage() ) ;
        }
        int currLevel = staxenv.getLevel();
        if (currLevel >=1)
        {
            ListIterator li = staxenv.getHandlerStackIterator(currLevel);
            while (li.hasPrevious())
            {
                Object ob = li.previous();
                if (ob instanceof AGAVEFeatureCallbackItf)
                {
                    ((AGAVEFeatureCallbackItf) ob).reportFeature(
                          new RangeLocation(region.getStart(), region.getEnd() ) );
                    return;
                }
            }

        }
  }
}