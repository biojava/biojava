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
package org.biojava.nbio.tornadovm;

import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.CifFileReader;

/** Test structures bundled as resources, so the tests need no network. */
public final class TestStructures {

	private TestStructures() {
	}

	public static Structure load(String id) throws Exception {
		try (InputStream in = new GZIPInputStream(
				TestStructures.class.getResourceAsStream("/pdb/" + id.toLowerCase() + ".cif.gz"))) {
			return new CifFileReader().getStructure(in);
		}
	}

	public static Atom[] caAtoms(String id, String chain) throws Exception {
		return StructureTools.getRepresentativeAtomArray(load(id).getPolyChainByPDB(chain));
	}

	/** Runs the given code with the GPU path forced on (or, without a TornadoVM runtime, the fallback). */
	public static <T> T forced(ThrowingSupplier<T> code) throws Exception {
		System.setProperty(TornadoSupport.PROPERTY, "force");
		try {
			return code.get();
		} finally {
			System.clearProperty(TornadoSupport.PROPERTY);
		}
	}

	public interface ThrowingSupplier<T> {
		T get() throws Exception;
	}
}
