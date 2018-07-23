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
package org.biojava.nbio.structure.io.mmtf;

import org.junit.Test;
import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.mmcif.ChemCompGroupFactory;
import org.biojava.nbio.structure.io.mmcif.DownloadChemCompProvider;

/**
 * Test bond finding in BioJava
 * @author Anthony Bradley
 *
 */
public class TestBondFinding {

	/**
	 * Test that the bonds we are finding is consistenty.
	 * @throws IOException
	 * @throws StructureException
	 */
	@Test
	public void testInterGroupBonds() throws IOException, StructureException {
		// Normal
		assertEquals(2236, getInterBonds("1QMZ"));
		// 	Disulphide
		assertEquals(956, getInterBonds("2QWO"));
		// Covalent ligand
		assertEquals(2294, getInterBonds("4QDV"));
		// DNA 
		assertEquals(22, getInterBonds("4XSN"));

	}

	/**
	 * Find all of the inter group bonds in a structure
	 * @param pdbId the pdb id of the structure to determine
	 * @return the number of inter group bonds (double counted) in a structure
	 * @throws IOException
	 * @throws StructureException
	 */
	public int getInterBonds(String pdbId) throws IOException, StructureException{
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		cache.setFetchBehavior(FetchBehavior.FETCH_FILES);
		FileParsingParameters params = cache.getFileParsingParams();
		params.setCreateAtomBonds(true);
		params.setAlignSeqRes(true);
		params.setParseBioAssembly(true);
		DownloadChemCompProvider dcc = new DownloadChemCompProvider();
		ChemCompGroupFactory.setChemCompProvider(dcc);
		dcc.checkDoFirstInstall();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
		int counter =0;
		// Now get the structure
		Structure newStruc = StructureIO.getStructure(pdbId);
		// Now loop through the atoms
		for(Chain c: newStruc.getChains()){
			for(Group g: c.getAtomGroups()){
				List<Atom> theseAtoms = g.getAtoms();
				for(Atom a: theseAtoms){
					List<Bond> theseBonds = a.getBonds();
					if(theseBonds != null){
						for(Bond b: a.getBonds()){
							Atom other = b.getOther(a);
							int indexOther = theseAtoms.indexOf(other);
							// Check if the index is within the group
							if(indexOther<0 || indexOther >= theseAtoms.size()){
								counter++;
							}
						}
					}
				}
			}
		}
		return counter;
	}
}
