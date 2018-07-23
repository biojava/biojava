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
package demo;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.cluster.SubunitClustererMethod;
import org.biojava.nbio.structure.cluster.SubunitClustererParameters;
import org.biojava.nbio.structure.gui.BiojavaJmol;
import org.biojava.nbio.structure.symmetry.axis.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;

import java.io.IOException;
import java.util.List;

/**
 * This demo shows how to display the {@link QuatSymmetryResults} of a
 * structure.
 * <p>
 * Examples: 4HHB, 4AQ5, 1LTI, 1STP, 4F88, 2W6E, 2LXC, 3OE7, 4INU, 4D8s, 4EAR,
 * 4IYQ, 3ZKR
 * <p>
 * Local symmetry: 2WPD (2 local symmetries), 4F88 (local C8), 1LTI (local C5),
 * 2W6E (local C3), 2LXC (local C2), 3OE7 (local C3)
 * <p>
 * Local Pseudosymmetry: 3ZDY, 3ZDX
 * <p>
 * Helical: 1B47
 * <p>
 * With internal symmetry: 4E3E, 1VYM
 * 
 * @author Peter Rose
 * @author Aleix Lafita
 * 
 */
public class DemoQuatSymmetryJmol {

	public static void main(String[] args) throws IOException,
			StructureException {

		String name = "2vml";

		// Download the biological assembly
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		Structure structure = cache.getStructure("BIO:" + name + ":1");

		QuatSymmetryParameters sp = new QuatSymmetryParameters();
		SubunitClustererParameters cp = new SubunitClustererParameters();
		cp.setClustererMethod(SubunitClustererMethod.SEQUENCE); // normal
		// cp.setClustererMethod(SubunitClustererMethod.STRUCTURE); // pseudo
		cp.setSequenceCoverageThreshold(0.9);

		// Calculate and display the global symmetry
		QuatSymmetryResults globalSymmetry = QuatSymmetryDetector
				.calcGlobalSymmetry(structure, sp, cp);
		showResults(structure, name, globalSymmetry);

		// Calculate and displaythe local symmetry
		List<QuatSymmetryResults> localSymmetry = QuatSymmetryDetector
				.calcLocalSymmetries(structure, sp, cp);

		for (QuatSymmetryResults result : localSymmetry)
			showResults(structure, name, result);

	}

	private static void showResults(Structure s, String name,
			QuatSymmetryResults results) {

		String title = name + ": " + results.getStoichiometry()
				+ ", " + results.getSymmetry();

		if (results.isPseudoStoichiometric())
			title += ", pseudosymmetric";

		if (results.isLocal())
			title += ", local";

		String script = "set defaultStructureDSSP true; set measurementUnits ANGSTROMS;  select all;  spacefill off; wireframe off; "
				+ "backbone off; cartoon on; color cartoon structure; color structure;  select ligand;wireframe 0.16;spacefill 0.5; "
				+ "color cpk ; select all; model 0;set antialiasDisplay true; autobond=false;save STATE state_1;";

		AxisAligner aligner = AxisAligner.getInstance(results);

		JmolSymmetryScriptGenerator scriptGenerator = JmolSymmetryScriptGeneratorPointGroup
				.getInstance(aligner, "g");

		script += scriptGenerator.getOrientationWithZoom(0);
		script += scriptGenerator.drawPolyhedron();
		script += scriptGenerator.drawAxes();
		script += scriptGenerator.colorBySymmetry();

		title += ", method: " + results.getMethod();

		script += "draw axes* on; draw poly* on;";

		BiojavaJmol jmol = new BiojavaJmol();
		jmol.setStructure(s);

		jmol.setTitle(title);
		jmol.evalString(script);

	}
}
