/**
 *  BioJava development code This code may be freely distributed and modified
 *  under the terms of the GNU Lesser General Public Licence. This should be
 *  distributed with the code. If you do not have a copy, see:
 *  http://www.gnu.org/copyleft/lesser.html Copyright for this code is held
 *  jointly by the individual authors. These should be listed in
 *
 *@author    doc comments. For more information on the BioJava project and its
 *      aims, or to join the biojava-l mailing list, visit the home page at:
 *      http://www.biojava.org/
 */

package org.biojava.bio.seq.io.game12;

import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.SAXException;

/**
 *  Handles the GAME &lt;dbxref&gt; element
 *
 * @author     David Huen
 * @since      1.2
 */
public class GAMEPropertyHandler
         extends StAXFeatureHandler {
    // <dbxref> is a container for external database references.
    // it is possible that non-unique <dbxref> occur.

    // temporary cache
    String propertyType = null;
    String propertyValue = null;

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_PROPERTY_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMEPropertyHandler(staxenv);
            }
        };


    /**
     *  Constructor for the GAMEDbxrefHandler object
     *
     *@param  staxenv   Description of the Parameter
     *@param  parentID  Description of the Parameter
     */
    GAMEPropertyHandler(StAXFeatureHandler staxenv) {
        // setup environment
        super(staxenv);

        // setup handlers
        // <type>
        super.addHandler(new ElementRecognizer.ByLocalName("type"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new TypeHandler();
                }
            }
                );
        // <value>
        super.addHandler(new ElementRecognizer.ByLocalName("value"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new ValueHandler();
                }
            }
                );
    }

    /**
     *  retrieves the &lt;type&gt; field.
     */
    private class TypeHandler extends StringElementHandlerBase {
        protected void setStringValue(String s) {
            propertyType = s.trim();
        }
    }

    /**
     *  retrieves the &lt;value&gt; field.
     */
    private class ValueHandler extends StringElementHandlerBase {
        protected void setStringValue(String s) {
            propertyValue = s.trim();
        }
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler) 
        throws SAXException 
    {
        // validate before going further
        if ((propertyType == null) || (propertyValue == null)) {
            return;
        }

        // set up field
        try {
            listener.addFeatureProperty(propertyType, propertyValue);
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("unexpected exception while adding <property> as a feature property.");
        }

    }
}
