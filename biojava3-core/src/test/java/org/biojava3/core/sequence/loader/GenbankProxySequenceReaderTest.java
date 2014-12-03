package org.biojava3.core.sequence.loader;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.features.FeatureInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Testing example for issue #834
 *
 * @author Jacek Grzebyta
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
    public void biojava3() throws Throwable {
        logger.info("run test for protein: {}", gi);
        GenbankProxySequenceReader<AminoAcidCompound> genbankReader
                = new GenbankProxySequenceReader<AminoAcidCompound>(System.getProperty("java.io.tmpdir"), 
                                                                    this.gi, 
                                                                    AminoAcidCompoundSet.getAminoAcidCompoundSet());

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
