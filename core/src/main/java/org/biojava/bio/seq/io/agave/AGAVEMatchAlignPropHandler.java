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
import org.biojava.utils.ChangeVetoException;
import org.xml.sax.SAXException;

/**
 * Deals with match_align
 * @author Hanning Ni    Doubletwist Inc
 */
public class AGAVEMatchAlignPropHandler
               extends StAXPropertyHandler
{
  // set up factory method
  public static final StAXHandlerFactory AGAVE_MATCH_ALIGN_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEMatchAlignPropHandler(staxenv);
    }
  };


  AGAVEMatchAlignPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("match_align", true);

  }

   public void characters(char[] ch, int start, int length)
        throws SAXException
  {
      try{
          staxenv.featureTemplate.annotation.setProperty("match_align", new String(ch) );
      }catch (ChangeVetoException cve) {
        throw new SAXException(" change veto exception ") ;
    }
  }

}
