package org.biojava.nbio.structure.symmetry.gui;

import java.awt.event.KeyEvent;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.MultipleAlignmentJmol;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;
import org.biojava.nbio.structure.align.util.RotationAxis;
import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryDetector;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryParameters;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.core.Subunits;
import org.biojava.nbio.structure.symmetry.internal.SymmetryAxes;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGenerator;
import org.biojava.nbio.structure.symmetry.jmolScript.JmolSymmetryScriptGeneratorPointGroup;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;

/**
 * Class that provides visualizations methods for symmetry
 * alignments. Call the display() method for the default 
 * visualization of symmetry.
 * 
 * @author Aleix Lafita
 * 
 */
public class SymmetryDisplay {

	/**
	 * Displays a multiple alignment of the symmetry subunits.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol displaySubunits(MultipleAlignment msa) 
			throws StructureException {

		MultipleAlignment subunits = SymmetryTools.toSubunitAlignment(msa);
		return MultipleAlignmentDisplay.display(subunits);
	}

	/**
	 * Displays a multiple alignment of the whole structure transformations
	 * colored by blocks, corresponding to the subunits.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol displayFull(MultipleAlignment msa) 
			throws StructureException {

		MultipleAlignment full = SymmetryTools.toFullAlignment(msa);

		MultipleAlignmentJmol jmol = MultipleAlignmentDisplay.display(full);
		jmol.setColorByBlocks(true);

		return jmol;
	}

	/**
	 * Displays a single structure in a cartoon representation with each
	 * symmetric subunit colored differently.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @param axes symmetry axes
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol display(MultipleAlignment msa,
			SymmetryAxes axes) throws StructureException {

		List<Atom[]> atoms = msa.getAtomArrays();

		MultipleAlignmentJmol jmol = new MultipleAlignmentJmol(msa, atoms);

		addSymmetryMenu(jmol, axes);

		//Show all the axes and point group symmetry
		if (axes!=null) jmol.evalString(printSymmetryAxes(msa, axes, false));
		jmol.evalString(printPointGroupAxes(msa));

		return jmol;
	}

	/**
	 * Displays a single structure in a cartoon representation with each
	 * symmetric subunit colored differently.
	 * 
	 * @param msa the symmetry multiple alignment obtained from CeSymm
	 * @throws StructureException
	 */
	public static MultipleAlignmentJmol display(MultipleAlignment msa)
			throws StructureException {
		return display(msa, null);
	}

	/**
	 * Adds a Symmetry menu to the Jmol display, so that further symmetry
	 * analysis can be triggered.
	 * 
	 * @param jmol parent jmol
	 * @param axes symmetry axes
	 */
	private static void addSymmetryMenu(MultipleAlignmentJmol jmol, 
			SymmetryAxes axes){

		JMenuBar menubar = jmol.getFrame().getJMenuBar();

		JMenu symm = new JMenu("Symmetry");
		symm.setMnemonic(KeyEvent.VK_S);

		SymmetryListener li = new SymmetryListener(jmol, axes);

		JMenuItem subunits = new JMenuItem("Subunit Superposition");
		subunits.addActionListener(li);
		symm.add(subunits);

		JMenuItem multiple = new JMenuItem("Multiple Structure Alignment");
		multiple.addActionListener(li);
		symm.add(multiple);

		JMenuItem pg = new JMenuItem("Point Group Symmetry");
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
	 * @param elementary only print elementary axes if true
	 * @return
	 */
	public static String printSymmetryAxes(MultipleAlignment msa, 
			SymmetryAxes axes, boolean elementary) {

		int id = 0;
		String script = "";
		Atom[] atoms = msa.getAtomArrays().get(0);

		List<Matrix4d> symmAxes = null;
		if (elementary){
			symmAxes = axes.getElementaryAxes();
		} else {
			symmAxes = axes.getSymmetryAxes();
		}

		for (Matrix4d axis : symmAxes) {
			RotationAxis rot = new RotationAxis(axis);
			script += rot.getJmolScript(atoms, id);
			id++;
		}
		return script;
	}

	/**
	 * Given a symmetry alignment, it draws the point group symmetry axes
	 * and the polyhedron box around the structure. 
	 * It uses the quaternary symmetry detection code, but tries to factor
	 * out the alignment and detection steps.
	 * 
	 * @param symm
	 * @return
	 */
	public static String printPointGroupAxes(MultipleAlignment symm){

		//Obtain the clusters of aligned Atoms and subunit variables
		MultipleAlignment subunits = SymmetryTools.toSubunitAlignment(symm);
		List<Atom[]> alignedCA = subunits.getAtomArrays();
		List<Integer> corePos = MultipleAlignmentTools.getCorePositions(
				subunits.getBlock(0));

		List<Point3d[]> caCoords = new ArrayList<Point3d[]>();
		List<Integer> folds = new ArrayList<Integer>();
		List<Boolean> pseudo = new ArrayList<Boolean>();
		List<String> chainIds = new ArrayList<String>();
		List<Integer> models = new ArrayList<Integer>();
		List<Double> seqIDmin = new ArrayList<Double>();
		List<Double> seqIDmax = new ArrayList<Double>();
		List<Integer> clusterIDs = new ArrayList<Integer>();
		int fold = 1;
		Character chain = 'A';

		for (int str=0; str<alignedCA.size(); str++){
			Atom[] array = alignedCA.get(str);
			List<Point3d> points = new ArrayList<Point3d>();
			List<Integer> alignedRes = 
					subunits.getBlock(0).getAlignRes().get(str);
			for (int pos=0; pos<alignedRes.size(); pos++){
				Integer residue = alignedRes.get(pos);
				if (residue == null) continue;
				else if (!corePos.contains(pos)) continue;
				Atom a = array[residue];
				points.add(new Point3d(a.getCoords()));
			}
			caCoords.add(points.toArray(new Point3d[points.size()]));
			if (alignedCA.size() % fold == 0){
				folds.add(fold); //the folds are the common denominators
			}
			fold++;
			pseudo.add(false);
			chainIds.add(chain+"");
			chain++;
			models.add(0);
			seqIDmax.add(1.0);
			seqIDmin.add(1.0);
			clusterIDs.add(0);
		}

		//Create directly the subunits, because we know the aligned CA
		Subunits globalSubunits = new Subunits(caCoords, clusterIDs, 
				pseudo, seqIDmin, seqIDmax, 
				folds, chainIds, models);

		//Quaternary Symmetry Detection
		QuatSymmetryParameters param = new QuatSymmetryParameters();
		param.setRmsdThreshold(symm.size() * 1.5);

		QuatSymmetryResults gSymmetry = 
				QuatSymmetryDetector.calcQuatSymmetry(globalSubunits, param);

		AxisAligner axes = AxisAligner.getInstance(gSymmetry);

		//Draw the axes as in the quaternary symmetry
		JmolSymmetryScriptGenerator scriptGenerator = 
				JmolSymmetryScriptGeneratorPointGroup.getInstance(axes, "g");

		String script = "save selection; set defaultStructureDSSP true; "
				+ "set measurementUnits ANGSTROMS;  select all;  "
				+ "spacefill off; wireframe off;"
				+ "set antialiasDisplay true; autobond=false; ";

		script += scriptGenerator.getOrientationWithZoom(0);
		script += "restore selection; ";
		script += scriptGenerator.drawPolyhedron();
		script += scriptGenerator.drawAxes();
		script += "draw axes* on; draw poly* on; ";

		return script;
	}

}
