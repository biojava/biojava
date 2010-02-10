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
 * <code>SSBindWUblastn2_0a19Test</code> tests object bindings for
 * Blast-like SAX events.
 *
 * @author Keith James
 * @since 1.2
 */
public class SSBindWUblastn2_0a19Test extends SSBindCase
{
    public SSBindWUblastn2_0a19Test(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
        super.setUp();

        setTopHitValues(12875d, "U51677",
                        1, 2575, StrandedFeature.POSITIVE,
                        1, 2575, StrandedFeature.POSITIVE);

        setBotHitValues(2123d, "D63874",
                        1, 2575, StrandedFeature.POSITIVE,
                        77, 974, StrandedFeature.POSITIVE);

        String resName = "org/biojava/bio/program/ssbind/wu_blastn_2.0a19.out.gz";
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

        assertEquals(4, result.getHits().size());
    }
}
