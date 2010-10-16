/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.DataSource;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenericFastaHeaderParserTest {

    public GenericFastaHeaderParserTest() {
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
     * GenBank                           gi|gi-number|gb|accession|locus
     * ENA Data Library                 gi|gi-number|emb|accession|locus
     * DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
     * NBRF PIR                          pir||entry
     * Protein Research Foundation       prf||name
     * SWISS-PROT                        sp|accession|name
     * Brookhaven Protein Data Bank (1)  pdb|entry|chain
     * Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
     * PDB EBI                           PDB:1ECY_A mol:protein length:142  ECOTIN
     * Patents                           pat|country|number
     * GenInfo Backbone Id               bbs|number
     * General database identifier       gnl|database|identifier
     * NCBI Reference Sequence           ref|accession|locus
     * Local Sequence identifier         lcl|identifier
     *
     * @author Scooter Willis <willishf at gmail dot com>
     */
    @Test
    public void testParseHeader() {
        System.out.println("parseHeader");
        String header = "";
        ProteinSequence sequence = new ProteinSequence("");
        GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound> instance = new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>();

        header = "gi|gi-number|gb|accession|locus";
        instance.parseHeader(header, sequence);
        assertEquals("accession", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.GENBANK);

        header = "gi|gi-number|emb|accession|locus";
        instance.parseHeader(header, sequence);
        assertEquals("accession", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.ENA);

        header = "gi|gi-number|dbj|accession|locus";
        instance.parseHeader(header, sequence);
        assertEquals("accession", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.DDBJ);

        header = "pir||entry";
        instance.parseHeader(header, sequence);
        assertEquals("entry", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.NBRF);

        header = "prf||name";
        instance.parseHeader(header, sequence);
        assertEquals("name", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.PRF);

        header = "sp|accession|name";
        instance.parseHeader(header, sequence);
        assertEquals("accession", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.UNIPROT);

        header = "pdb|entry|chain";
        instance.parseHeader(header, sequence);
        assertEquals("entry:chain", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.PDB1);

        header = "entry:chain|PDBID|CHAIN|SEQUENCE";
        instance.parseHeader(header, sequence);
        assertEquals("entry:chain", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.PDB2);
        header = "PDB:1ECY_A mol:protein length:142  ECOTIN";
        instance.parseHeader(header, sequence);
        assertEquals("1ECY_A", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.PDBe);

        header = "pat|country|number";
        instance.parseHeader(header, sequence);
        assertEquals("number", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.PATENTS);

        header = "bbs|number";
        instance.parseHeader(header, sequence);
        assertEquals("number", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.GENINFO);

        header = "gnl|database|identifier";
        instance.parseHeader(header, sequence);
        assertEquals("identifier", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.GENERAL);

        header = "ref|accession|locus";

        instance.parseHeader(header, sequence);
        assertEquals("accession", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.NCBI);

        header = "lcl|identifier";
        instance.parseHeader(header, sequence);
        assertEquals("identifier", sequence.getAccession().getID());
        assertEquals(sequence.getAccession().getDataSource(), DataSource.LOCAL);
        // TODO review the generated test code and remove the default call to fail.
        //fail("The test case is a prototype.");
    }
}
