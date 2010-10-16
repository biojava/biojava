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

import org.xml.sax.SAXException;

/**
 * handle AGAVE xref
 * @author Hanning Ni ,      Doubletwist Inc
 */
public class AGAVEXrefPropHandler  extends StAXPropertyHandler
     implements AGAVEDbIdCallbackItf, AGAVEDbIdPropCallbackItf
{
   // set up factory method
  public static final StAXHandlerFactory AGAVE_XREF_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEXrefPropHandler(staxenv);
    }
  };

  private AGAVEXref xref ;

  AGAVEXrefPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("xref", true);
    xref = new AGAVEXref() ;
    super.addHandler(new ElementRecognizer.ByLocalName("db_id"),
         AGAVEDbIdPropHandler.AGAVE_DBID_PROP_HANDLER_FACTORY);
    super.addHandler(new ElementRecognizer.ByLocalName("xref_property"),
         AGAVEXrefPropPropHandler.AGAVE_XREF_PROP_PROP_HANDLER_FACTORY);

  }

  public void addDbId(AGAVEDbId db_id)
  {
      xref.setDbId( db_id) ;
  }

  public void addProperty( AGAVEProperty prop)
  {
     xref.addProp( prop ) ;
  }

   public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler handler)
                throws SAXException
  {
        int currLevel = staxenv.getLevel();
        if (currLevel >=1)
        {
            ListIterator li = staxenv.getHandlerStackIterator(currLevel);
            while (li.hasPrevious())
            {
                Object ob = li.previous();
                if (ob instanceof AGAVEXrefCallbackItf)
                {
                    ((AGAVEXrefCallbackItf) ob).addXref( xref );
                    return;
                }
            }

        }
  }
}
