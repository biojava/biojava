/*
 *                  BioJava development code
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
 * Created on Jun 8, 2007
 *
 */
package org.biojava.bio.structure;

import java.io.IOException;
import java.io.InputStream;

import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.io.PDBFileParser;

import junit.framework.TestCase;

public class StructurePairAlignerTest extends TestCase {

    Structure structure1;
    Structure structure2;

    protected void setUp()
    {
        InputStream inStream = this.getClass().getResourceAsStream("/5pti.pdb");
        assertNotNull(inStream);
        InputStream inStream2 = this.getClass().getResourceAsStream("/1tap.pdb");
        assertNotNull(inStream2);

        PDBFileParser pdbpars = new PDBFileParser();
        try {
            structure1 = pdbpars.parsePDBFile(inStream) ;
            structure2 = pdbpars.parsePDBFile(inStream2);
        } catch (IOException e) {
            e.printStackTrace();
        }

        assertNotNull(structure1);
        assertNotNull(structure2);
        assertEquals("structure does not contain one chain ", 1 ,structure1.size());

    }


    public void testAlignStructureStructure() {
       StructurePairAligner aligner = new StructurePairAligner();

       boolean allFine = true;
       String msg = "";
       try {
           aligner.align(structure1,structure2);

           AlternativeAlignment[] aligs = aligner.getAlignments();
           assertEquals("the number of obtained alternative alignments is not correct",20, aligs.length);
           AlternativeAlignment a = aligs[0];
           
           assertNotNull(a);
           
           assertEquals("the expected nr of eq. residues is not correct.",47,a.getEqr());
           
           // they are v. close, but not identical
           assertTrue(a.getRmsd() < 4);
           assertTrue(a.getRmsd() > 3);
           assertTrue(a.getPercId() > 9);
           assertTrue(a.getScore() > 140);

       } catch (Exception e){
           msg = e.getMessage();
           allFine = false;
           e.printStackTrace();
       }
       assertTrue(allFine);
       assertEquals("an error occured","" ,msg);



    }



}
