/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.features.FeatureInterface;
import org.biojava.nbio.core.sequence.features.PublicationReference;
import org.biojava.nbio.core.sequence.features.PublicationReferenceAuthor;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.junit.*;

import java.io.InputStream;
import java.util.LinkedHashMap;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class GenBankReferenceTest {

    public GenBankReferenceTest() {
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


    @Test
    public void testDirectSubmission() throws Exception {

        InputStream inStream = this.getClass().getResourceAsStream("/SCU49845.gb");

        GenbankReader<DNASequence, NucleotideCompound> GenbankDNA =
                new GenbankReader<DNASequence, NucleotideCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
                        new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
                );

        LinkedHashMap<String, DNASequence> dnaSequences = GenbankDNA.process();
        DNASequence result = dnaSequences.get("U49845");
        List<PublicationReference> publicationReferences = result.getPublicationReference();
        Assert.assertEquals(2, publicationReferences.size());
        PublicationReference publicationReference1 = publicationReferences.get(0);
        Assert.assertEquals("8846915", publicationReference1.getId());
        Assert.assertEquals("Genes Dev. 10 (7), 777-793 (1996)", publicationReference1.getJournal());
        Assert.assertEquals(PublicationReference.ReferenceType.PUBMED, publicationReference1.getReferenceType());
        Assert.assertEquals("Selection of axial growth sites in yeast requires Axl2p, a novelplasma membrane glycoprotein", publicationReference1.getTitle());

        List<PublicationReferenceAuthor> authors = publicationReference1.getAuthors();
        Assert.assertEquals(4, authors.size());

        PublicationReferenceAuthor author1 = authors.get(0);
        Assert.assertEquals("T.", author1.getFirstName());
        Assert.assertEquals("Roemer", author1.getLastName());
        Assert.assertEquals("T. Roemer", author1.getFullName());

        PublicationReferenceAuthor author2 = authors.get(1);
        Assert.assertEquals("K.", author2.getFirstName());
        Assert.assertEquals("Madden", author2.getLastName());
        Assert.assertEquals("K. Madden", author2.getFullName());

        PublicationReferenceAuthor author3 = authors.get(2);
        Assert.assertEquals("J.", author3.getFirstName());
        Assert.assertEquals("Chang", author3.getLastName());
        Assert.assertEquals("J. Chang", author3.getFullName());

        PublicationReferenceAuthor author4 = authors.get(3);
        Assert.assertEquals("M.", author4.getFirstName());
        Assert.assertEquals("Snyder", author4.getLastName());
        Assert.assertEquals("M. Snyder", author4.getFullName());

        PublicationReference publicationReference2 = publicationReferences.get(1);
        //Direct submission have no id
        Assert.assertEquals(null, publicationReference2.getId());
        Assert.assertEquals("Submitted (22-FEB-1996) Biology, Yale University, New Haven, CT06520, USA", publicationReference2.getJournal());
        Assert.assertEquals(PublicationReference.ReferenceType.DIRECT_SUBMISSION, publicationReference2.getReferenceType());
        Assert.assertEquals("Direct Submission", publicationReference2.getTitle());

        inStream.close();
    }

    @Test
    public void testPatent() throws Exception {

        InputStream inStream = this.getClass().getResourceAsStream("/E01172.gb");

        GenbankReader<DNASequence, NucleotideCompound> GenbankDNA =
                new GenbankReader<DNASequence, NucleotideCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
                        new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
                );


        LinkedHashMap<String, DNASequence> dnaSequences = GenbankDNA.process();

        DNASequence result = dnaSequences.get("E01172");
        List<PublicationReference> publicationReferences = result.getPublicationReference();
        Assert.assertEquals(1, publicationReferences.size());
        PublicationReference publicationReference1 = publicationReferences.get(0);
        Assert.assertEquals(null, publicationReference1.getId());
        Assert.assertEquals("Patent: JP 1987099398-A 1 08-MAY-1987;F HOFFMANN LA ROCHE & CO AG", publicationReference1.getJournal());
        Assert.assertEquals(PublicationReference.ReferenceType.PATENT, publicationReference1.getReferenceType());
        Assert.assertEquals("NOVEL POLYPEPTIDE", publicationReference1.getTitle());
        List<PublicationReferenceAuthor> authors = publicationReference1.getAuthors();

        Assert.assertEquals(4, authors.size());

        inStream.close();
    }

    @Test
    public void testMore() throws Exception {

        InputStream inStream = this.getClass().getResourceAsStream("/HE608876.gb");

        GenbankReader<DNASequence, NucleotideCompound> GenbankDNA =
                new GenbankReader<DNASequence, NucleotideCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
                        new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
                );

        LinkedHashMap<String, DNASequence> dnaSequences = GenbankDNA.process();

        assertEquals(1, dnaSequences.size());
        DNASequence result = dnaSequences.values().iterator().next();

        inStream.close();
    }

    @Test
    public void testKeywords() throws Exception {

        InputStream inStream = this.getClass().getResourceAsStream("/NM_000266.gb");

        GenbankReader<DNASequence, NucleotideCompound> GenbankDNA =
                new GenbankReader<DNASequence, NucleotideCompound>(
                        inStream,
                        new GenericGenbankHeaderParser<DNASequence, NucleotideCompound>(),
                        new DNASequenceCreator(DNACompoundSet.getDNACompoundSet())
                );

        LinkedHashMap<String, DNASequence> dnaSequences = GenbankDNA.process();

        assertEquals(1, dnaSequences.size());
        DNASequence result = dnaSequences.values().iterator().next();
        List<FeatureInterface<AbstractSequence<NucleotideCompound>, NucleotideCompound>> features = result.getFeatures();

        inStream.close();
    }

}