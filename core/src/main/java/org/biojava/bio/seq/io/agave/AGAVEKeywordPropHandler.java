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
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.SAXException;

/**
 * Deals with AGAVE keywords
 *
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEKeywordPropHandler
               extends StAXPropertyHandler
{
  // set up factory method
  public static final StAXHandlerFactory AGAVE_KEYWORD_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEKeywordPropHandler(staxenv);
    }
  };


  AGAVEKeywordPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("keyword", true);
  }

   public void characters(char[] ch, int start, int length)
        throws SAXException
  {
      try{
         String kw = (String) UtilHelper.getProperty(staxenv.featureTemplate.annotation,"keyword");
         if( kw != null )
             kw = kw + "<keyword>" +  new String( ch ) + "</keyword>" + "\n" ;
         else
         {
             kw = "<keyword>" +  new String( ch ) + "</keyword>" + "\n" ;
         }
         staxenv.featureTemplate.annotation.setProperty("keyword", kw);
      }catch (ChangeVetoException cve) {
        throw new SAXException(" change veto exception ") ;
    }
  }

}
