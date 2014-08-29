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

/**
 * Testing example for issue #834
 *
 * @author Jacek Grzebyta
 * @see InfoTask
 */
@RunWith(Parameterized.class)
public class GenbankProxySequenceReaderTest {

    private String gi;

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
    public void biojava3Approach1() throws Throwable {

        GenbankProxySequenceReader<AminoAcidCompound> genbankReader
                = new GenbankProxySequenceReader<AminoAcidCompound>("/tmp", this.gi, AminoAcidCompoundSet.getAminoAcidCompoundSet());

        ProteinSequence seq = new ProteinSequence(genbankReader, AminoAcidCompoundSet.getAminoAcidCompoundSet());

        Assert.assertNotNull("protein sequence is null", seq);
        genbankReader.getHeaderParser().parseHeader(genbankReader.getHeader(), seq);
        genbankReader.getFeatureParser().parseFeatures(seq);
        
        Assert.assertTrue(seq.getDescription() != null);

        Assert.assertFalse(seq.getFeaturesKeyWord().getKeyWords().isEmpty());
        Assert.assertFalse(seq.getFeaturesByType("organism").get(0).getSource().isEmpty());
        
        System.err.println(String.format("taxonomy id: %s", seq.getTaxonomy()));
        Assert.assertNotNull(seq.getTaxonomy().getID());
        Assert.assertNotNull(seq.getSequenceAsString());
        
        
        List<FeatureInterface<AbstractSequence<AminoAcidCompound>, AminoAcidCompound>> codedBy = seq.getFeaturesByType("coded_by");

        if (!codedBy.isEmpty()) {
            // get parent DNA
            Assert.assertNotNull(seq.getParentSequence().getSequenceAsString() != null);
        }
    }
}
