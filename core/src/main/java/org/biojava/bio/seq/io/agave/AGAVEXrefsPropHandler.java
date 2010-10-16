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

import org.biojava.utils.ChangeVetoException;
import org.xml.sax.SAXException;

/**
 * Deals with database crossreferences (xrefs)
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEXrefsPropHandler
               extends StAXPropertyHandler implements AGAVEDbIdCallbackItf, AGAVEXrefCallbackItf
{
  // set up factory method
  public static final StAXHandlerFactory AGAVE_XREFS_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEXrefsPropHandler(staxenv);
    }
  };
  private AGAVEXrefs xrefs ;

   AGAVEXrefsPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("xrefs", true);
    xrefs = new AGAVEXrefs() ;
    super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);
    super.addHandler(new ElementRecognizer.ByLocalName("xref"),
         AGAVEXrefPropHandler.AGAVE_XREF_PROP_HANDLER_FACTORY);

  }

  public void addDbId(AGAVEDbId db_id)
  {
      xrefs.addDbId( db_id) ;
  }

  public void addXref( AGAVEXref xref)
  {
      xrefs.addXref( xref ) ;
  }

   public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
       try{
          List set = (List)UtilHelper.getProperty(staxenv.featureTemplate.annotation,"xrefs") ;
          if( set == null )
          {
              set = new ArrayList(1);
              set.add( xrefs );
              staxenv.featureTemplate.annotation.setProperty( "xrefs", set) ;
          }
          else
          {
              set.add( xrefs );
          }
       }catch(ChangeVetoException e){
         throw new SAXException( e.getMessage() ) ;
       }
  }
}
