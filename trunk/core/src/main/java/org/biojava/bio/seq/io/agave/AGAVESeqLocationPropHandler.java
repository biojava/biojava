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

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.RangeLocation;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * seq_location
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVESeqLocationPropHandler
 extends StAXPropertyHandler
{

   // set up factory metho
  public static final StAXHandlerFactory AGAVE_SEQ_LOCATION_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVESeqLocationPropHandler(staxenv);
    }
  };
   private int start , end ;
   private String strand = "false" ;

   AGAVESeqLocationPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("seq_location", true);
   }

  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
      start =new Integer( attrs.getValue( "least_start" ) ) .intValue() ;
      end = new Integer( attrs.getValue( "greatest_end" ) ) .intValue() ;
      strand = attrs.getValue("is_on_complement") ;
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
		// THOMASD fixed strand == null crash

                Object ob = li.previous();
                if (ob instanceof AGAVEFeatureCallbackItf)
                {
                    ((AGAVEFeatureCallbackItf) ob).reportFeature(  new RangeLocation(start, end) );
                    if( "true".equalsIgnoreCase(strand) )
                        ((AGAVEFeatureCallbackItf) ob).reportStrand( StrandedFeature.NEGATIVE );
                   else if( "false".equalsIgnoreCase(strand) )
                         ((AGAVEFeatureCallbackItf) ob).reportStrand( StrandedFeature.POSITIVE );
                   else
                       ((AGAVEFeatureCallbackItf) ob).reportStrand( StrandedFeature.UNKNOWN );
                    return;
                }
            }

        }
   }
}
