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

import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava.bio.structure.io.PDBFileParser;

import junit.framework.TestCase;

public class StructureToolsTest extends TestCase {

    Structure structure;

    protected void setUp()
    {
        InputStream inStream = this.getClass().getResourceAsStream("/5pti.pdb");
        assertNotNull(inStream);


        PDBFileParser pdbpars = new PDBFileParser();
        FileParsingParameters params = new FileParsingParameters();
        params.setAlignSeqRes(false);
        pdbpars.setFileParsingParameters(params);
        
        try {
            structure = pdbpars.parsePDBFile(inStream) ;
        } catch (IOException e) {
            e.printStackTrace();
        }

        assertNotNull(structure);

        assertEquals("structure does not contain one chain ", 1 ,structure.size());
    }


    public void testGetCAAtoms(){
        Atom[] cas = StructureTools.getAtomCAArray(structure);
        assertEquals("did not find the expected number of Atoms (58), but got " + cas.length,58,cas.length);
    }

    public void testGetNrAtoms(){
        int length = StructureTools.getNrAtoms(structure);
        assertEquals("did not find the expected number of Atoms (1070), but got " + length,1070,length);


    }



}
