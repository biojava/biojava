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

import java.io.BufferedInputStream;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.biojava.bio.program.sax.BlastLikeSAXParser;
import org.biojava.bio.search.SeqSimilaritySearchResult;
import org.biojava.bio.seq.StrandedFeature;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;

/**
 * <code>SSBindWUtblastn2_0a19Test</code> tests object bindings for
 * Blast-like SAX events.
 *
 * @author Keith James
 * @since 1.2
 */
public class SSBindWUtblastn2_0a19Test extends SSBindCase
{
    public SSBindWUtblastn2_0a19Test(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        super.setUp();

        setTopHitValues(1144d, "X12597",
                        1, 214, null,
                        56, 697, StrandedFeature.POSITIVE);

        setBotHitValues(126d, "AF025664",
                        47, 214, null,
                        304, 866, StrandedFeature.POSITIVE);

        String resName = "org/biojava/bio/program/ssbind/wu_tblastn_2.0a19.out.gz";
        InputStream resStream = getClass().getClassLoader().getResourceAsStream(
                resName);
        assert resStream != null
                : "Resource " + resName + " could not be located";
        searchStream =
            new GZIPInputStream(new BufferedInputStream(resStream));

        // XMLReader -> (SAX events) -> adapter -> builder -> objects
        XMLReader reader = (XMLReader) new BlastLikeSAXParser();

        reader.setContentHandler(adapter);
        reader.parse(new InputSource(searchStream));
    }

    public void testBlastResultHitCount()
    {
        SeqSimilaritySearchResult result =
            (SeqSimilaritySearchResult) searchResults.get(0);

        assertEquals(250, result.getHits().size());
    }
}
