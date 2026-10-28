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
package org.biojava.nbio.tornadovm.bench;

import java.io.File;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.io.CifFileReader;

/** Loads benchmark structures from {@code data/pdb/<id>.cif.gz} when present, else through {@link StructureIO}. */
final class Structures {

	private Structures() {
	}

	static Structure load(String id) throws Exception {
		File f = new File(System.getProperty("bench.pdbdir", "data/pdb"), id.toLowerCase() + ".cif.gz");
		if (f.isFile()) {
			return new CifFileReader().getStructure(f);
		}
		return StructureIO.getStructure(id);
	}
}
