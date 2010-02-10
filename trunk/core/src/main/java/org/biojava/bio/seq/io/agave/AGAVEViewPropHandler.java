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
 * Moves view attributes into annotation properties.
 * <p>
 * start --> view_start
 * <p>
 * end --> view_end
 *
 * @author Hanning Ni    Doubletwist Inc
 * @author Keith James
 * @author Matthew Pocock
 */
public class AGAVEViewPropHandler
               extends StAXPropertyHandler
{
  // set up factory method
  public static final StAXHandlerFactory AGAVE_VIEW_PROP_HANDLER_FACTORY
    = new StAXHandlerFactory() {
    public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
      return new AGAVEViewPropHandler(staxenv);
    }
  };

  AGAVEViewPropHandler(StAXFeatureHandler staxenv) {
    // execute superclass method to setup environment
    super(staxenv);
    setHandlerCharacteristics("view", true);

  }
  public void startElementHandler(
                String nsURI,
                String localName,
                String qName,
                Attributes attrs)
         throws SAXException
  {
      try{
        staxenv.featureTemplate.annotation.setProperty("view_start",  attrs.getValue("start") ) ;
        staxenv.featureTemplate.annotation.setProperty("view_length", attrs.getValue("length") ) ;
      }catch(Exception e){
        throw new SAXException( e.getMessage() ) ;
      }
  }
}

