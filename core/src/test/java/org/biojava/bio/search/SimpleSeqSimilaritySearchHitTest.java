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

package org.biojava.bio.search;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.SimpleAlignment;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.SimpleSymbolList;

/**
 * <code>SimpleSeqSimilaritySearchHitTest</code> tests the behaviour
 * of <code>SimpleSeqSimilaritySearchHit</code>.
 *
 * @author Keith James
 */
public class SimpleSeqSimilaritySearchHitTest extends TestCase
{
    private SeqSimilaritySearchHit h1;
    private SeqSimilaritySearchHit h2;

    private List subHits1;
    private List subHits2;

    private Alignment al1;
    private Alignment al2;

    private double            score = 100.0d;
    private double           eValue = 1e-10d;
    private double           pValue = 1e-10d;
    private int          queryStart = 1;
    private int            queryEnd = 10;
    private Strand   querySeqStrand = StrandedFeature.POSITIVE;
    private int        subjectStart = 2;
    private int          subjectEnd = 8;
    private Strand subjectSeqStrand = StrandedFeature.POSITIVE;

    private String        subjectID = "subjectID";
    private String   querySeqTokens = "TRYPASNDEF";
    private String subjectSeqTokens = "-RYPASND--";

    public SimpleSeqSimilaritySearchHitTest(String name)
    {
        super(name);
    }

    protected void setUp() throws Exception
    {
	SymbolTokenization tp = ProteinTools.getAlphabet().getTokenization("token");

        Map labelMap1 = new HashMap();
        labelMap1.put(SeqSimilaritySearchSubHit.QUERY_LABEL,
                      new SimpleSymbolList(tp, querySeqTokens));
        labelMap1.put(subjectID,
                      new SimpleSymbolList(tp, subjectSeqTokens));

        al1 = new SimpleAlignment(labelMap1);

        Map labelMap2 = new HashMap();
        labelMap2.put(SeqSimilaritySearchSubHit.QUERY_LABEL,
                      new SimpleSymbolList(tp, querySeqTokens));
        labelMap2.put(subjectID,
                      new SimpleSymbolList(tp, subjectSeqTokens));

        al2 = new SimpleAlignment(labelMap2);

        subHits1 = new ArrayList();
        subHits1.add(new SimpleSeqSimilaritySearchSubHit(score,
                                                         eValue,
                                                         pValue,
                                                         queryStart,
                                                         queryEnd,
                                                         querySeqStrand,
                                                         subjectStart,
                                                         subjectEnd,
                                                         subjectSeqStrand,
                                                         al1,
                                                         Annotation.EMPTY_ANNOTATION));

        subHits2 = new ArrayList();
        subHits2.add(new SimpleSeqSimilaritySearchSubHit(score,
                                                         eValue,
                                                         pValue,
                                                         queryStart,
                                                         queryEnd,
                                                         querySeqStrand,
                                                         subjectStart,
                                                         subjectEnd,
                                                         subjectSeqStrand,
                                                         al2,
                                                         Annotation.EMPTY_ANNOTATION));

         h1 = new SimpleSeqSimilaritySearchHit(score,
                                               eValue,
                                               pValue,
                                               queryStart,
                                               queryEnd,
                                               querySeqStrand,
                                               subjectStart,
                                               subjectEnd,
                                               subjectSeqStrand,
                                               subjectID,
                                               Annotation.EMPTY_ANNOTATION,
                                               subHits1);

         h2 = new SimpleSeqSimilaritySearchHit(score,
                                               eValue,
                                               pValue,
                                               queryStart,
                                               queryEnd,
                                               querySeqStrand,
                                               subjectStart,
                                               subjectEnd,
                                               subjectSeqStrand,
                                               subjectID,
                                               Annotation.EMPTY_ANNOTATION,
                                               subHits2);
    }

    public void testEquals()
    {
        assertEquals(h1, h1);
        assertEquals(h2, h2);
        assertEquals(h1, h2);
        assertEquals(h2, h1);
    }

    public void testScores()
    {
        assertEquals(h1.getScore(),  100.0d, 0.0d);
        assertEquals(h1.getEValue(), 1e-10d, 0.0d);
        assertEquals(h1.getPValue(), 1e-10d, 0.0d);
    }

    public void testQuery()
    {
        assertEquals(h1.getQueryStart(), 1);
        assertEquals(h1.getQueryEnd(),  10);
        assertEquals(h1.getQueryStrand(), StrandedFeature.POSITIVE);
    }

    public void testSubject()
    {
        assertEquals(h1.getSubjectStart(), 2);
        assertEquals(h1.getSubjectEnd(),   8);
        assertEquals(h1.getSubjectStrand(), StrandedFeature.POSITIVE);
    }

    public void testID()
    {
        assertEquals(h1.getSubjectID(), "subjectID");
    }

    public void testAnnotation()
    {
        assertEquals(h1.getAnnotation(), Annotation.EMPTY_ANNOTATION);
    }
}
