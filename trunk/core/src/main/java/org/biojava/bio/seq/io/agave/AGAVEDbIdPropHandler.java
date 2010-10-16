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
 * Deals with database crossreferences
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEDbIdPropHandler
               extends StAXPropertyHandler
{

   public static final StAXHandlerFactory AGAVE_DBID_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEDbIdPropHandler(staxenv);
    }
   };

   private AGAVEDbId db_id ;
   AGAVEDbIdPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("db_id", true);
  }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
      db_id = new AGAVEDbId() ;
      db_id.setId( attrs.getValue( "id" ) );
      db_id.setVersion( attrs.getValue( "version" ) );
      db_id.setDbCode( attrs.getValue( "db_code" ) ) ;
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
              if( ob instanceof AGAVEDbIdCallbackItf )
              {
                  ( (AGAVEDbIdCallbackItf) ob ).addDbId( db_id ) ;
                 //  return ;
              }
           }
       }
       try{
       //if no suitable handler found, save it to annotation
        staxenv.featureTemplate.annotation.setProperty( "db_id", db_id) ;
       }catch(Exception e){
          throw new SAXException( e.getMessage() ) ;
      }
  }


}
