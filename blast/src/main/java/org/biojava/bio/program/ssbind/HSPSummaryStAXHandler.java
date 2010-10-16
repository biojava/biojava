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

import org.biojava.bio.search.SearchContentHandler;
import org.biojava.utils.stax.StAXContentHandler;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;

/**
 * <code>HSPSummaryStAXHandler</code> handles the HSPSummary element
 * of BioJava BlastLike XML.
 *
 * @author Keith James
 * @since 1.3
 */
public class HSPSummaryStAXHandler extends SeqSimilarityStAXHandler
{
    public static final StAXHandlerFactory HSPSUMMARY_HANDLER_FACTORY =
        new StAXHandlerFactory()
        {
            public StAXContentHandler getHandler(SeqSimilarityStAXAdapter ssContext)
            {
                return new HSPSummaryStAXHandler(ssContext);
            }
        };

    HSPSummaryStAXHandler(SeqSimilarityStAXAdapter ssContext)
    {
        super(ssContext);
    }

    protected void handleStartElement(String     nsURI,
                                      String     localName,
                                      String     qName,
                                      Attributes attrs)
        throws SAXException
    {
        SearchContentHandler sch = ssContext.getSearchContentHandler();

        // score, expectValue, numberOfIdentities, alignmentSize and
        // percentageIdentity always present in valid XML
        sch.addSubHitProperty("score", attrs.getValue("score"));
        sch.addSubHitProperty("expectValue", attrs.getValue("expectValue"));
        sch.addSubHitProperty("numberOfIdentities",
                              attrs.getValue("numberOfIdentities"));
        sch.addSubHitProperty("alignmentSize",
                              attrs.getValue("alignmentSize"));
        sch.addSubHitProperty("percentageIdentity",
                              attrs.getValue("percentageIdentity"));

        String attr;
        attr = attrs.getValue("bitScore");
        if (attr != null)
            sch.addSubHitProperty("bitScore", attr);

        attr = attrs.getValue("pValue");
        if (attr != null)
            sch.addSubHitProperty("pValue", attr);

        attr = attrs.getValue("numberOfIdentities");
        if (attr != null)
            sch.addSubHitProperty("numberOfIdentities", attr);

        attr = attrs.getValue("numberOfPositives");
        if (attr != null)
            sch.addSubHitProperty("numberOfPositives", attr);

        attr = attrs.getValue("querySequenceType");
        if (attr != null)
            sch.addSubHitProperty("querySequenceType", attr);

        attr = attrs.getValue("hitSequenceType");
        if (attr != null)
            sch.addSubHitProperty("subjectSequenceType", attr);

        // These are not explicitly set by BLASTP
        attr = attrs.getValue("queryStrand");
        if (attr != null)
            sch.addSubHitProperty("queryStrand", attr);

        attr = attrs.getValue("hitStrand");
        if (attr != null)
            sch.addSubHitProperty("subjectStrand", attr);

        attr = attrs.getValue("queryFrame");
        if (attr != null)
            sch.addSubHitProperty("queryFrame", attr);

        attr = attrs.getValue("hitFrame");
        if (attr != null)
            sch.addSubHitProperty("subjectFrame", attr);
    }
}
