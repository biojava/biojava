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
package org.biojava.nbio.structure.symmetry.gui;

import java.awt.event.KeyEvent;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentJmolDisplay;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;
import org.jmol.util.Logger;

/**
 * Class that provides visualizations methods for symmetry alignments. Call the
 * display() method for the default visualization of symmetry.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryDisplay {

	/**
	 * Displays a multiple alignment of the symmetry repeats.
	 *
	 * * @param symm CeSymmResult
	 *
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol displayRepeats(CeSymmResult symm)
			throws StructureException {

		MultipleAlignment repeats = SymmetryTools.toRepeatsAlignment(symm);
		MultipleAlignmentJmol jmol = MultipleAlignmentJmolDisplay.display(repeats);
		jmol.setTitle(getSymmTitle(symm));
		return jmol;
	}

	/**
	 * Displays a multiple alignment of the whole structure transformations
	 * colored by blocks, corresponding to the symmetric protodomains.
	 *
	 * @param symm
	 *            CeSymmResult
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol displayFull(CeSymmResult symm)
			throws StructureException {

		MultipleAlignment full = SymmetryTools.toFullAlignment(symm);

		MultipleAlignmentJmol jmol = MultipleAlignmentJmolDisplay.display(full);
		jmol.setColorByBlocks(true);
		jmol.setTitle(getSymmTitle(symm));

		return jmol;
	}

	/**
	 * Displays a single structure in a cartoon representation with each
	 * symmetric repeat colored differently.
	 *
	 * @param msa
	 *            the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static AbstractAlignmentJmol display(CeSymmResult symmResult)
			throws StructureException {

		if (symmResult.isSignificant() && symmResult.isRefined()) {
			// Show the structure colored by repeat
			MultipleAlignment msa = symmResult.getMultipleAlignment();
			List<Atom[]> atoms = msa.getAtomArrays();
			MultipleAlignmentJmol jmol = new MultipleAlignmentJmol(msa, atoms);
			jmol.setTitle(jmol.getStructure().getPDBHeader().getTitle());
			addSymmetryMenu(jmol, symmResult);
			jmol.evalString(printSymmetryGroup(symmResult));
			jmol.evalString(printSymmetryAxes(symmResult, false));
			jmol.setTitle(getSymmTitle(symmResult));
			return jmol;
		} else {
			// Show the optimal self-alignment
			Logger.info("Showing optimal self-alignment");
			Atom[] cloned = StructureTools
					.cloneAtomArray(symmResult.getAtoms());
			AbstractAlignmentJmol jmol = StructureAlignmentDisplay.display(
					symmResult.getSelfAlignment(), symmResult.getAtoms(),
					cloned);
			RotationAxis axis = new RotationAxis(symmResult.getSelfAlignment());
			jmol.evalString(axis.getJmolScript(symmResult.getAtoms()));
			jmol.setTitle(getSymmTitle(symmResult));
			return jmol;
		}
	}

	/**
	 * Adds a Symmetry menu to the Jmol display, so that further symmetry
	 * analysis can be triggered.
	 *
	 * @param jmol
	 *            parent jmol
	 * @param symmResult
	 *            CeSymmResult
	 */
	private static void addSymmetryMenu(MultipleAlignmentJmol jmol,
			CeSymmResult symmResult) {

		JMenuBar menubar = jmol.getFrame().getJMenuBar();

		JMenu symm = new JMenu("Symmetry");
		symm.setMnemonic(KeyEvent.VK_S);

		SymmetryListener li = new SymmetryListener(jmol, symmResult);

		JMenuItem repeats = new JMenuItem("Repeats Superposition");
		repeats.addActionListener(li);
		symm.add(repeats);

		JMenuItem multiple = new JMenuItem("Multiple Structure Alignment");
		multiple.addActionListener(li);
		symm.add(multiple);

		JMenuItem self = new JMenuItem("Optimal Self Alignment");
		self.addActionListener(li);
		symm.add(self);

		JMenuItem pg = new JMenuItem("Show Symmetry Group");
		pg.addActionListener(li);
		symm.add(pg);

		JMenuItem ax = new JMenuItem("Show Symmetry Axes");
		ax.addActionListener(li);
		symm.add(ax);

		JMenuItem news = new JMenuItem("New Symmetry Analysis");
		news.addActionListener(li);
		symm.add(news);

		menubar.add(symm, 3);
		jmol.getFrame().pack();
	}

	/**
	 * Generates a String that displays the symmetry axes of a structure.
	 *
	 * @param msa
	 * @param axes
	 * @param elementary
	 *            only print elementary axes if true
	 * @return
	 */
	public static String printSymmetryAxes(CeSymmResult symm, boolean elementary) {

		int id = 0;
		String script = "";
		Atom[] atoms = symm.getAtoms();

		List<Matrix4d> symmAxes = null;
		if (elementary) {
			symmAxes = symm.getAxes().getElementaryAxes();
		} else {
			symmAxes = symm.getAxes().getSymmetryAxes();
		}

		for (Matrix4d axis : symmAxes) {
			RotationAxis rot = new RotationAxis(axis);
			script += rot.getJmolScript(atoms, id);
			id++;
		}
		return script;
	}

	/**
	 * Given a symmetry alignment, it draws the symmetry group axes and the
	 * polyhedron box around the structure. It uses the quaternary symmetry
	 * detection code, but tries to factor out the alignment and detection
	 * steps.
	 *
	 * @param symm
	 *            CeSymmResult
	 * @return
	 * @throws StructureException
	 */
	public static String printSymmetryGroup(CeSymmResult symm) throws StructureException {

		QuatSymmetryResults gSymmetry = SymmetryTools.getQuaternarySymmetry(symm);

		AxisAligner axes = AxisAligner.getInstance(gSymmetry);

		// Draw the axes as in the quaternary symmetry
		JmolSymmetryScriptGenerator scriptGenerator = JmolSymmetryScriptGeneratorPointGroup
				.getInstance(axes, "g");

		String script = "save selection; set measurementUnits ANGSTROMS;"
				+ "select all; set antialiasDisplay true; autobond=false; ";

		script += scriptGenerator.getInstantaneousOrientation(0);
		script += "restore selection; ";
		script += scriptGenerator.drawPolyhedron();
		script += scriptGenerator.drawAxes();
		script += "draw axes* on; draw poly* on; ";

		return script;
	}

	/**
	 * Create a symmetry title for a display frame (Jmol, alignment, etc). The
	 * title contains information about the algorithm, structure id and
	 * parameters used.
	 *
	 * @param result
	 * @return title String
	 */
	public static String getSymmTitle(CeSymmResult result) {

		StringBuffer buff = new StringBuffer();

		// Add algorithm name and version
		buff.append(result.getMultipleAlignment().getEnsemble()
				.getAlgorithmName());
		buff.append(" V");
		buff.append(result.getMultipleAlignment().getEnsemble().getVersion());
		buff.append(": ");

		// Add the result summary string
		buff.append(result.toString());
		return buff.toString();
	}

}
