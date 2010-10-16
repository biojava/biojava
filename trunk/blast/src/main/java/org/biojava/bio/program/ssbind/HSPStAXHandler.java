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
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <code>HSPStAXHandler</code> handles the HSP element of BioJava
 * BlastLike XML.
 *
 * @author Keith James
 * @since 1.3
 */
public class HSPStAXHandler extends SeqSimilarityStAXHandler
{
    public static final StAXHandlerFactory HSP_HANDLER_FACTORY =
        new StAXHandlerFactory()
        {
            public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
            {
                return new HSPStAXHandler(ssContext);
            }
        };

    /**
     * Creates a new instance which sends callbacks to the specified
     * <code>SeqSimilarityStAXAdapter</code>.
     *
     * @param ssContext a <code>SeqSimilarityStAXAdapter</code>.
     */
    HSPStAXHandler(SeqSimilarityStAXAdapter ssContext)
    {
        super(ssContext);
        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "HSPSummary"),
                   HSPSummaryStAXHandler.HSPSUMMARY_HANDLER_FACTORY);

        addHandler(new ElementRecognizer.ByNSName(SeqSimilarityStAXAdapter.NAMESPACE,
                                                  "BlastLikeAlignment"),
                   AlignmentStAXHandler.ALIGNMENT_HANDLER_FACTORY);
    }

    protected void handleStartElement(String     nsURI,
                                      String     localName,
                                      String     qName,
                                      Attributes attrs)
        throws SAXException
    {
        ssContext.getSearchContentHandler().startSubHit();
    }

    protected void handleEndElement(String     nsURI,
                                    String     localName,
                                    String     qName)
        throws SAXException
    {
        ssContext.getSearchContentHandler().endSubHit();
    }
}
