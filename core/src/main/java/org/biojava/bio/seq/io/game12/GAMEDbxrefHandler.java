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
public class GAMEDbxrefHandler
         extends StAXFeatureHandler {
    // <dbxref> is a container for external database references.
    // it is possible that non-unique <dbxref> occur.

    // temporary cache
    String db_xref_db;
    String db_xref_id;

    // set up factory method
    /**
     *  Description of the Field
     */
    public final static StAXHandlerFactory GAME_DBXREF_HANDLER_FACTORY
             =
        new StAXHandlerFactory() {
            public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                return new GAMEDbxrefHandler(staxenv);
            }
        };


    /**
     *  Constructor for the GAMEDbxrefHandler object
     *
     *@param  staxenv   Description of the Parameter
     *@param  parentID  Description of the Parameter
     */
    GAMEDbxrefHandler(StAXFeatureHandler staxenv) {
        // setup environment
        super(staxenv);

        // setup handlers
        // <xref_db> : external database name
        super.addHandler(new ElementRecognizer.ByLocalName("xref_db"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new XrefDbHandler();
                }
            }
                );
        // <db_xref_id> : external database id
        super.addHandler(new ElementRecognizer.ByLocalName("db_xref_id"),
            new StAXHandlerFactory() {
                public StAXContentHandler getHandler(StAXFeatureHandler staxenv) {
                    return new DbXrefIdHandler();
                }
            }
                );
    }

    /**
     *  retrieves the &lt;xref_db&gt; field.
     */
    private class XrefDbHandler extends StringElementHandlerBase {
        protected void setStringValue(String s) {
            db_xref_db = s.trim();
        }
    }

    /**
     *  retrieves the &lt;db_xref_id&gt; field.
     */
    private class DbXrefIdHandler extends StringElementHandlerBase {
        protected void setStringValue(String s) {
            db_xref_id = s.trim();
        }
    }

    // does not have its own returnData() as it does not expect
    // to have any returned to it.

    public void endElementHandler(
            String nsURI,
            String localName,
            String qName,
            StAXContentHandler contentHandler) 
        throws SAXException 
    {
        // validate before going further
        if ((db_xref_db == null) || (db_xref_id == null)) {
            return;
        }
    }
}
