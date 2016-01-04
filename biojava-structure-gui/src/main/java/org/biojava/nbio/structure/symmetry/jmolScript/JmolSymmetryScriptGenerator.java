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
package org.biojava.nbio.structure.symmetry.jmolScript;

import org.biojava.nbio.structure.symmetry.core.AxisAligner;
import org.biojava.nbio.structure.symmetry.core.HelixAxisAligner;
import org.biojava.nbio.structure.symmetry.core.RotationAxisAligner;

import javax.vecmath.Color4f;
import javax.vecmath.Matrix4d;
import javax.vecmath.Tuple3d;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

public abstract class JmolSymmetryScriptGenerator {

	/**
	 * Returns an instance of a JmolSymmetryScriptGenerator, based on the symmetry of a structure (factory method)
	 * @param axisAligner
	 * @param rotationGroup
	 * @return instance of JmolSymmetryScriptGenerator
	 */
	public static JmolSymmetryScriptGenerator getInstance(AxisAligner axisAligner, String name) {
		String symmetry = axisAligner.getSymmetry();
		
		if (symmetry.equals("C1")) {
			return new JmolSymmetryScriptGeneratorC1((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.startsWith("C")) {
			return new JmolSymmetryScriptGeneratorCn((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.startsWith("D")) {
			return new JmolSymmetryScriptGeneratorDn((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("T")) {
			return new JmolSymmetryScriptGeneratorT((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("O")) {
			return new JmolSymmetryScriptGeneratorO((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("I")) {
			return new JmolSymmetryScriptGeneratorI((RotationAxisAligner)axisAligner, name);
		} else if (symmetry.equals("H")) {
			return new JmolSymmetryScriptGeneratorH((HelixAxisAligner)axisAligner, name);
		}
		
		return null;
	}
	/**
	 * Returns the Jmol zoom to fit polyhedron and symmetry axes. This zoom
	 * level should be used so that the polyhedron and symmetry axes are not cutoff.
	 * @return
	 */
	abstract public int getZoom();

	/**
	 * Returns a Jmol script to set the default orientation for a structure
	 * @return Jmol script
	 */
	public abstract String getDefaultOrientation();

	/**
	 * Returns the number of orientations available for this structure
	 * @return number of orientations
	 */
	public abstract int getOrientationCount();

	/**
	 * Returns a Jmol script that sets a specific orientation
	 * @param index orientation index
	 * @return Jmol script
	 */
	public abstract String getOrientation(int index);
	
	/**
	 * Returns a Jmol script that sets a specific orientation instantaneously
	 * @param index orientation index
	 * @return Jmol script
	 */
	public String getInstantaneousOrientation(int index){
		String s = getOrientation(index);
		return s.replaceAll("moveto 4", "moveto 0");
	}

	/**
	 * Returns a Jmol script that sets a specific orientation and zoom
	 * to draw either axes or polyhedron
	 * @param index orientation index
	 * @return Jmol script
	 */
	public abstract String getOrientationWithZoom(int index);

	/**
	 * Returns the name of a specific orientation
	 * @param index orientation index
	 * @return name of orientation
	 */
	public abstract String getOrientationName(int index);

	/**
	 * Returns transformation matrix to orient structure
	 * @return transformation matrix
	 */
	public abstract Matrix4d getTransformation();
	
	/** Sets a default Jmol script used for coloring. This method is
	 * used in local symmetry cases to color those subunits that are
	 * not related by symmetry.
	 * @param colorScript
	 */	
	public abstract void setDefaultColoring(String colorScript);
	
	/**
	 * Sets the type of bioassembly to be colored. If set to true,
	 * it will generate a Jmol script for a bioassembly generated
	 * by Jmol on the fly. If set to false, it will generate Jmol script for
	 * a bioassembly file read by Jmol.
	 */
	public abstract void setOnTheFly(boolean onTheFly);
	
	/**
	 * Returns a Jmol script that draws an invisible polyhedron around a structure.
	 * Use showPolyhedron() and hidePolyhedron() to toggle visibility.
	 * @return Jmol script
	 */
	public abstract String drawPolyhedron();

	public abstract String hidePolyhedron();

	public abstract String showPolyhedron();

	/**
	 * Returns a Jmol script that draws symmetry or inertia axes for a structure.
	 * Use showAxes() and hideAxes() to toggle visibility.
	 * @return Jmol script
	 */
	public abstract String drawAxes();

	/**
	 * Returns a Jmol script to hide axes
	 * @return Jmol script
	 */
	public abstract String hideAxes();

	/**
	 * Returns a Jmol script to show axes
	 * @return Jmol script
	 */
	public abstract String showAxes();

	/**
	 * Returns a Jmol script that displays a symmetry polyhedron and symmetry axes
	 * and then loop through different orientations
	 * @return Jmol script
	 */
	public abstract String playOrientations();

	/**
	 * Returns a Jmol script that colors the subunits of a structure by different colors
	 * @return
	 */
	public abstract String colorBySubunit();

	/**
	 * Returns a Jmol script that colors subunits by their sequence cluster ids.
	 * @return Jmol script
	 */
	public abstract String colorBySequenceCluster();

	/**
	 * Returns a Jmol script that colors subunits to highlight the symmetry within a structure
	 * @return Jmol script
	 */
	public abstract String colorBySymmetry();
	
	protected static String getJmolColorScript(Map<Color4f, List<String>> map) {
		StringBuilder s = new StringBuilder();
		s.append("color cartoons none;");
		for (Entry<Color4f, List<String>> entry: map.entrySet()) {
			s.append("color{");
			List<String> ids = entry.getValue();
			for (int i = 0; i < ids.size(); i++) {
				s.append(":");
				s.append(ids.get(i));
				if (i < ids.size() -1 ) {
			    s.append("|");
				}
			}
			s.append("}");	
			s.append(getJmolColor(entry.getKey()));
			s.append(";");		
		}
		return s.toString();
	}
	
	protected static String getJmolColor(Color4f color) {
        String hex = Integer.toHexString((color.get().getRGB() & 0xffffff) | 0x1000000).substring(1);
        return " [x" + hex + "]";
	}
	
	protected static String getJmolPoint(Tuple3d point) {
		StringBuilder s = new StringBuilder();
		s.append("{");
		s.append(fDot2(point.x));
		s.append(",");
		s.append(fDot2(point.y));
		s.append(",");
		s.append(fDot2(point.z));
		s.append("}");
		return s.toString();
	}
	
	protected static String f1Dot2(float number) {
		return String.format("%1.2f", number);
	}
	
	protected static String fDot2(double number) {
		return String.format("%.2f", number);
	}
	
	/**
	 * Returns a lower precision floating point number for Jmol
	 * @param f
	 * @return
	 */
	protected static float jMolFloat(double f) {
		return (float)f;
	}
	
	protected static String getJmolLigandScript() {
		return "select ligand;wireframe 0.16;spacefill 23%;color cpk;";
	}

}