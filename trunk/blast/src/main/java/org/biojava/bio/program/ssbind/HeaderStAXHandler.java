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

package org.biojava.bio.program.ssbind;

import org.biojava.bio.program.xff.ElementRecognizer;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <code>HeaderStAXHandler</code> handles the Header element of
 * BioJava BlastLike XML.
 *
 * @author Keith James
 * @author Matthew Pocock
 * @since 1.3
 */
public class HeaderStAXHandler extends SeqSimilarityStAXHandler
{
    public static final StAXHandlerFactory HEADER_HANDLER_FACTORY =
        new StAXHandlerFactory()
        {
            public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
            {
                return new HeaderStAXHandler(ssContext);
            }
        };

    /**
     * Creates a new instance which sends callbacks to the specified
     * <code>SeqSimilarityStAXAdapter</code>.
     *
     * @param ssContext a <code>SeqSimilarityStAXAdapter</code>.
     */
    HeaderStAXHandler(SeqSimilarityStAXAdapter ssContext)
    {
        super(ssContext);

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "QueryId"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new QueryIDStAXHandler();
                       }
                   });

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "QueryDescription"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new QueryDescriptionStAXHandler();
                       }
                   });

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "DatabaseId"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new DatabaseIDStAXHandler();
                       }
                   });
    }

    /**
     * <code>QueryIDStAXHandler</code> handles the query sequence ID.
     */
    private class QueryIDStAXHandler extends StAXContentHandlerBase
    {
        public void startElement(String            uri,
                                 String            localName,
                                 String            qName,
                                 Attributes        attr,
                                 DelegationManager dm)
        throws SAXException
        {
            ssContext.getSearchContentHandler().setQueryID(attr.getValue("id"));
            if (attr.getValue("queryLength") != null)
            {
            	ssContext.getSearchContentHandler().addSearchProperty("queryLength",
            			attr.getValue("queryLength"));
            }
        }
    }

    /**
     * <code>QueryDescriptionStAXHandler</code> handles the hit
     * description.
     */
    private class QueryDescriptionStAXHandler extends StringElementHandlerBase
    {
        protected void setStringValue(String s)
        {
            ssContext.getSearchContentHandler().addSearchProperty("queryDescription", s);
        }
    }

    /**
     * <code>DatabaseIDStAXHandler</code> handles the database ID.
     */
    private class DatabaseIDStAXHandler extends StAXContentHandlerBase
    {
        public void startElement(String            uri,
                                 String            localName,
                                 String            qName,
                                 Attributes        attr,
                                 DelegationManager dm)
        throws SAXException
        {
            ssContext.getSearchContentHandler().setDatabaseID(attr.getValue("id"));
        }
    }
}
