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

package org.biojava.bio.seq.io.game12;

import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.game.ElementRecognizer;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.SAXException;

/**
 * Handles the GAME &lt;aspect&gt; element
 *
 * @author David Huen
 * @since 1.8
 */
public class GAMEAspectHandler extends StAXFeatureHandler {
    // this class has <function> and <dbxref> child elements.
    private String aspectDbxrefDb = null;
    private String aspectDbxrefId = null;
    private String function = null;

    // set up factory method
    public static final StAXHandlerFactory GAME_ASPECT_HANDLER_FACTORY 
      = new StAXHandlerFactory() {
      public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
        return new GAMEAspectHandler(staxenv);
      }
    };

    private class DbxrefHandler extends GAMEDbxrefHandler
    {
        private DbxrefHandler(StAXFeatureHandler staxenv)
        {
            super(staxenv);
        }

        public void endElementHandler(
                String nsURI,
                String localName,
                String qName,
                StAXContentHandler contentHandler)
            throws SAXException
        {
            // validate before going further
            super.endElementHandler(nsURI, localName, qName, contentHandler);

            // stash away result until processing complete
            aspectDbxrefDb = db_xref_db;
            aspectDbxrefId = db_xref_id;
        }
    }

    private class FunctionHandler extends StringElementHandlerBase {
        /**
         *  Sets the stringValue attribute of the NameHandler object
         *
         *@param  s  The new stringValue value
         */
        protected void setStringValue(String s) {
            function = s.trim();
        }
    }

    GAMEAspectHandler(StAXFeatureHandler staxenv) {
        // execute superclass method to setup environment
        super(staxenv);

        // setup handlers
        // <function>
        super.addHandler(new ElementRecognizer.ByLocalName("function"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new FunctionHandler();
                }
            }
        );
        // <dbxref>
        super.addHandler(new ElementRecognizer.ByLocalName("dbxref"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new DbxrefHandler(staxenv);
                }
            }
        );
    }

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler)
        throws SAXException
    {
        if ((aspectDbxrefDb == null) || (aspectDbxrefId == null) || (function == null)) return;

        try {
            listener.addFeatureProperty("aspect:" + function, "dbxref:" + aspectDbxrefDb + "//" + aspectDbxrefId);
        }
        catch (ParseException pe) {
            pe.printStackTrace();
            throw new SAXException("unexpected exception while adding <aspect> as a feature property.");
        }
    }
}

