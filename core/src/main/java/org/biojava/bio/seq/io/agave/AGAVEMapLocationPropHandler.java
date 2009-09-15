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
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEMapLocationPropHandler  extends StAXPropertyHandler
{

   // set up factory method
  public static final StAXHandlerFactory AGAVE_MAP_LOCATION_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEMapLocationPropHandler(staxenv);
    }
  };
  private AGAVEMapLocation ml ;

  AGAVEMapLocationPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("map_location", true);
    ml = new AGAVEMapLocation() ;
    super.addHandler(new ElementRecognizer.ByLocalName("map_position"),
         AGAVEMapPositionPropHandler.AGAVE_MAP_POSITION_PROP_HANDLER_FACTORY);
  }

  public void addMapPosition(AGAVEMapPosition pos)
  {
      ml.addPosition( pos ) ;
  }
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
     ml.setMapType( attrs.getValue("map_type")  ) ;
     ml.setChromosome( attrs.getValue("chromsome")  ) ;
     ml.setUnits( attrs.getValue("units")  ) ;
     ml.setSource( attrs.getValue("source")  ) ;
     ml.setSubSeqStart( attrs.getValue("sebseq_start")  ) ;
     ml.setOrientation( attrs.getValue("orientation")  ) ;
  }



   public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
      try{
         Object ob =UtilHelper.getProperty(staxenv.featureTemplate.annotation,"map_location");
         if( ob != null )
             ((List)ob).add( ml ) ;
         else
         {
             List mls = new ArrayList(1) ;
             mls.add( ml ) ;
             staxenv.featureTemplate.annotation.setProperty("map_location", mls);
         }
      }catch(Exception e){
         throw new SAXException(e.getMessage() ) ;
      }
  }
}

