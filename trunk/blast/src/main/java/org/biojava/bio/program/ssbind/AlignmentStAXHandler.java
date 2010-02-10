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
import org.biojava.utils.stax.StringElementHandlerBase;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <code>AlignmentStAXHandler</code> handles the BlastLikeAlignment
 * element of BioJava BlastLike XML.
 *
 * @author Keith James
 * @since 1.3
 */
public class AlignmentStAXHandler extends SeqSimilarityStAXHandler
{
    public static final StAXHandlerFactory ALIGNMENT_HANDLER_FACTORY =
        new StAXHandlerFactory()
        {
            public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
            {
                return new AlignmentStAXHandler(ssContext);
            }
        };

    /**
     * Creates a new instance which sends callbacks to the specified
     * <code>SeqSimilarityStAXAdapter</code>.
     *
     * @param ssContext a <code>SeqSimilarityStAXAdapter</code>.
     */
    AlignmentStAXHandler(SeqSimilarityStAXAdapter ssContext)
    {
        super(ssContext);

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "QuerySequence"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new QuerySequenceStAXHandler();
                       }
                   });

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "HitSequence"),
                   new StAXHandlerFactory()
                   {
                       public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
                       {
                           return new HitSequenceStAXHandler();
                       }
                   });
    }

    /**
     * <code>QuerySequenceStAXHandler</code> handles the query
     * sequence.
     */
    private class QuerySequenceStAXHandler extends StringElementHandlerBase
    {
        private SearchContentHandler sch;

        public void startElement(String            nsURI,
                                 String            localName,
                                 String            qName,
                                 Attributes        attrs,
                                 DelegationManager dm)
            throws SAXException
        {
            super.startElement(nsURI, localName, qName, attrs, dm);

            sch = ssContext.getSearchContentHandler();
            sch.addSubHitProperty("querySequenceStart",
                                  attrs.getValue("startPosition"));
            sch.addSubHitProperty("querySequenceEnd",
                                  attrs.getValue("stopPosition"));
        }

        protected void setStringValue(String s) throws SAXException
        {
            sch = ssContext.getSearchContentHandler();
            sch.addSubHitProperty("querySequence", s);
        }
    }

    /**
     * <code>HitSequenceStAXHandler</code> handles the hit sequence.
     */
    private class HitSequenceStAXHandler extends StringElementHandlerBase
    {
        private SearchContentHandler sch;

        public void startElement(String            nsURI,
                                 String            localName,
                                 String            qName,
                                 Attributes        attrs,
                                 DelegationManager dm)
            throws SAXException
        {
            super.startElement(nsURI, localName, qName, attrs, dm);

            sch = ssContext.getSearchContentHandler();
            sch.addSubHitProperty("subjectSequenceStart",
                                  attrs.getValue("startPosition"));
            sch.addSubHitProperty("subjectSequenceEnd",
                                  attrs.getValue("stopPosition"));
        }

        protected void setStringValue(String s) throws SAXException
        {
            sch = ssContext.getSearchContentHandler();
            sch.addSubHitProperty("subjectSequence", s);
        }
    }
}
