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
import org.biojava.bio.search.SearchContentHandler;
import org.biojava.utils.stax.DelegationManager;
import org.biojava.utils.stax.StAXContentHandler;
import org.biojava.utils.stax.StAXContentHandlerBase;
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <code>AlignmentStAXHandler</code> handles the Hit element of
 * BioJava BlastLike XML.
 *
 * @author Keith James
 * @since 1.3
 */
public class HitStAXHandler extends SeqSimilarityStAXHandler
{
    public static final StAXHandlerFactory HIT_HANDLER_FACTORY =
        new StAXHandlerFactory()
        {
            public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
            {
                return new HitStAXHandler(ssContext);
            }
        };

    /**
     * Creates a new instance which sends callbacks to the specified
     * <code>SeqSimilarityStAXAdapter</code>.
     *
     * @param ssContext a <code>SeqSimilarityStAXAdapter</code>.
     */
    HitStAXHandler(SeqSimilarityStAXAdapter ssContext)
    {
        super(ssContext);
        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "HitId"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new HitIDStAXHandler();
                       }
                   });

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
                                                  "HitDescription"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new HitDescriptionStAXHandler();
                       }
                   });

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "HSPCollection"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new HSPCollectionStAXHandler();
                       }
                   });
    }

    protected void handleStartElement(String     nsURI,
                                      String     localName,
                                      String     qName,
                                      Attributes attrs)
        throws SAXException
    {
        SearchContentHandler sch = ssContext.getSearchContentHandler();

        sch.startHit();
        if (attrs.getValue("sequenceLength") != null)
        {
            sch.addHitProperty("subjectSequenceLength",
                               attrs.getValue("sequenceLength"));
        }
    }

    protected void handleEndElement(String     nsURI,
                                    String     localName,
                                    String     qName)
        throws SAXException
    {
        ssContext.getSearchContentHandler().endHit();
    }

    /**
     * <code>HitIDStAXHandler</code> handles the hit ID.
     */
    private class HitIDStAXHandler extends StAXContentHandlerBase
    {
        public void startElement(String            uri,
                                 String            localName,
                                 String            qName,
                                 Attributes        attr,
                                 DelegationManager dm)
        throws SAXException
        {
            ssContext.getSearchContentHandler().addHitProperty("subjectId", attr.getValue("id"));
        }
    }

    /**
     * <code>QueryIDStAXHandler</code> handles the query ID.
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
        }
    }

    /**
     * <code>HitDescriptionStAXHandler</code> handles the hit
     * description.
     */
    private class HitDescriptionStAXHandler extends StringElementHandlerBase
    {
        protected void setStringValue(String s)
        {
            ssContext.getSearchContentHandler().addHitProperty("subjectDescription", s);
        }
    }

    /**
     * <code>HSPCollectionStAXHandler</code> handles the HSPCollection
     * element.
     */
    private class HSPCollectionStAXHandler extends StAXContentHandlerBase
    {
        public void startElement(String            uri,
                                 String            localName,
                                 String            qName,
                                 Attributes        attr,
                                 DelegationManager dm)
        throws SAXException
        {
            if (localName.equals("HSP"))
                dm.delegate(HSPStAXHandler.HSP_HANDLER_FACTORY.getHandler(ssContext));
        }
    }
}
