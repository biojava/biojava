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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.sequence.io;

import static org.junit.Assert.assertNotNull;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.GenBankQualifierMap;
import org.biojava.nbio.core.sequence.features.Qualifier;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 * @author Jacek Grzebyta
 */
public class GenbankReaderTest {

    private final static Logger logger = LoggerFactory.getLogger(GenbankReaderTest.class);

    public GenbankReaderTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of process method, of class GenbankReader.
     */
    @Test
    public void testProcess() throws Exception {

        logger.info("process protein");
        InputStream inStream = this.getClass().getResourceAsStream("/BondFeature.gb");
        assertNotNull(inStream);

        GenbankReader<ProteinSequence, AminoAcidCompound> GenbankProtein
                = new GenbankReader<ProteinSequence, AminoAcidCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<ProteinSequence, AminoAcidCompound>(),
                        new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
                );
        @SuppressWarnings("unused")
        LinkedHashMap<String, ProteinSequence> proteinSequences = GenbankProtein.process();
        inStream.close();

        logger.info("process DNA");
        inStream = this.getClass().getResourceAsStream("/NM_000266.gb");
        assertNotNull(inStream);

        GenbankReader<DNASequence, NucleotideCompound> GenbankDNA
                = new GenbankReader<DNASequence, NucleotideCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
                        new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
                );
        @SuppressWarnings("unused")
        LinkedHashMap<String, DNASequence> dnaSequences = GenbankDNA.process();
        inStream.close();
    }


    @Test
    public void CDStest() throws Exception {
        logger.info("CDS Test");

        InputStream inStream = this.getClass().getResourceAsStream("/BondFeature.gb");
        assertNotNull(inStream);
        
        GenbankReader<ProteinSequence, AminoAcidCompound> GenbankProtein
                = new GenbankReader<ProteinSequence, AminoAcidCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<ProteinSequence, AminoAcidCompound>(),
                        new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet())
                );
        LinkedHashMap<String, ProteinSequence> proteinSequences = GenbankProtein.process();
        inStream.close();


        Assert.assertTrue(proteinSequences.size() == 1);
        logger.info("protein sequences: {}", proteinSequences);
        
        ProteinSequence protein = new ArrayList<ProteinSequence>(proteinSequences.values()).get(0);
        
        FeatureInterface<AbstractSequence<AminoAcidCompound>, AminoAcidCompound> cdsFeature = protein.getFeaturesByType("CDS").get(0);
        String codedBy = cdsFeature.getQualifierMap().get("coded_by").getFirstValue();
        GenBankQualifierMap quals = cdsFeature.getQualifierMap();
        Qualifier dbref = quals.get("db_xref");
        String[] dbrefValues=dbref.getValues();

        Assert.assertNotNull(codedBy);
        Assert.assertTrue(!codedBy.isEmpty());
        Assert.assertEquals(codedBy, "NM_000266.2:503..904");
        Assert.assertEquals(5, dbrefValues.length);

    }

}
