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
package org.biojava.nbio.core.sequence.loader;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * Testing example for issue #834
 *
 * @author Jacek Grzebyta
 * @author Paolo Pavan
 * @see InfoTask
 */
@RunWith(Parameterized.class)
public class GenbankProxySequenceReaderTest {

    private String gi;
    private final static Logger logger = LoggerFactory.getLogger(GenbankProxySequenceReaderTest.class);

    public GenbankProxySequenceReaderTest(String gi) {
        this.gi = gi;
    }

    @Parameterized.Parameters
    public static Collection<String[]> getExamples() {
        String[][] out = new String[][]{
            {"399235158"},
            {"7525057"},
            {"34567"},
            {"379015144"},
            {"381353147"},
            {"381353148"},
            {"152970917"},
            {"7856885"},
            {"381353149"},
            {"254839678"}
        };

        return Arrays.asList(out);
    }

    @Test
    public void biojava3() throws IOException, InterruptedException, CompoundNotFoundException  { 
        logger.info("run test for protein: {}", gi);
        GenbankProxySequenceReader<AminoAcidCompound> genbankReader
                = new GenbankProxySequenceReader<AminoAcidCompound>(System.getProperty("java.io.tmpdir"), 
                                                                    this.gi, 
                                                                    AminoAcidCompoundSet.getAminoAcidCompoundSet());

        // why only tests on protein sequences?
        ProteinSequence seq = new ProteinSequence(genbankReader, AminoAcidCompoundSet.getAminoAcidCompoundSet());

        Assert.assertNotNull("protein sequence is null", seq);
        genbankReader.getHeaderParser().parseHeader(genbankReader.getHeader(), seq);
        
        Assert.assertTrue(seq.getDescription() != null);

        Assert.assertFalse(seq.getFeaturesKeyWord().getKeyWords().isEmpty());
        Assert.assertFalse(seq.getFeaturesByType("source").get(0).getSource().isEmpty());
        
        logger.info("taxonomy id: {}", seq.getTaxonomy().getID());
        Assert.assertNotNull(seq.getTaxonomy().getID());
        Assert.assertNotNull(seq.getSequenceAsString());
        
        
        List<FeatureInterface<AbstractSequence<AminoAcidCompound>, AminoAcidCompound>> codedBy = seq.getFeaturesByType("coded_by");

        if (!codedBy.isEmpty()) {
            // get parent DNA
            Assert.assertNotNull(seq.getParentSequence().getSequenceAsString() != null);
        }
    }
}
