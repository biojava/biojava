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
/**
 * 
 */
package org.biojava.nbio.structure.symmetry.jmolScript;

import org.biojava.nbio.structure.symmetry.core.Rotation;
import org.biojava.nbio.structure.symmetry.core.RotationAxisAligner;
import org.biojava.nbio.structure.symmetry.core.RotationGroup;
import org.biojava.nbio.structure.symmetry.core.Subunits;
import org.biojava.nbio.structure.symmetry.geometry.Polyhedron;
import org.jcolorbrewer.ColorBrewer;

import javax.vecmath.*;

import java.awt.*;
import java.util.*;
import java.util.List;

/**
 * @author Peter
 *
 */
public abstract class JmolSymmetryScriptGeneratorPointGroup extends JmolSymmetryScriptGenerator {
	private static String N_FOLD_AXIS_COLOR = "red";
	private static String TWO_FOLD_AXIS_COLOR = "deepskyblue";
	private static String THREE_FOLD_AXIS_COLOR = "lime";
	private static double AXIS_SCALE_FACTOR = 1.2;
	
	private RotationAxisAligner rotationAxisAligner = null;
	private RotationGroup rotationGroup = null;
	private Polyhedron polyhedron = null;
	private String name = "";
	private String defaultColoring = "";
	private boolean onTheFly = true;
	
	public JmolSymmetryScriptGeneratorPointGroup(RotationAxisAligner rotationAxisAligner, String name) {
		this.rotationAxisAligner = rotationAxisAligner;
		this.rotationGroup = this.rotationAxisAligner.getRotationGroup();
		this.name = name;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getZoom()
	 */
	@Override
	abstract public int getZoom();
	
	public void setOnTheFly(boolean onTheFly) {
		this.onTheFly = onTheFly;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getDefaultOrientation()
	 */
	@Override
	public String getDefaultOrientation() {	
		StringBuilder s = new StringBuilder();
		s.append(setCentroid());
		
		// calculate  orientation
		Quat4d q = new Quat4d();
		q.set(polyhedron.getViewMatrix(0));
		q.normalize();
		Quat4d r = new Quat4d();
		r.set(rotationAxisAligner.getRotationMatrix());
		r.normalize();
		q.mul(r);
		q.normalize();
		
		// set orientation
		s.append("moveto 0 quaternion{");
		s.append(jMolFloat(q.x));
		s.append(",");
		s.append(jMolFloat(q.y));
		s.append(",");
		s.append(jMolFloat(q.z));
		s.append(",");
		s.append(jMolFloat(q.w));
		s.append("};");
		return s.toString();
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getOrientationCount()
	 */
	@Override
	public int getOrientationCount() {
		return polyhedron.getViewCount();
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getOrientation(int)
	 */
	@Override
	public String getOrientation(int index) {	
		StringBuilder s = new StringBuilder();
		s.append(setCentroid());
		
		// calculate  orientation
		Quat4d q = new Quat4d();
		q.set(polyhedron.getViewMatrix(index));
		q.normalize();
		Quat4d r = new Quat4d();
		r.set(rotationAxisAligner.getRotationMatrix());
		r.normalize();
		q.mul(r);
		q.normalize();
		
		// set orientation
		s.append("moveto 4 quaternion{");
		s.append(jMolFloat(q.x));
		s.append(",");
		s.append(jMolFloat(q.y));
		s.append(",");
		s.append(jMolFloat(q.z));
		s.append(",");
		s.append(jMolFloat(q.w));
		s.append("}");
		s.append(";");
		return s.toString();
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getOrientationWithZoom(int)
	 */
	@Override
	public String getOrientationWithZoom(int index) {
		StringBuilder s = new StringBuilder();
		s.append(getOrientation(index));
		s.insert(s.length()-1, getZoom());
		return s.toString();
		
	}
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getOrientationName(int)
	 */
	@Override
	public String getOrientationName(int index) {	
	    return polyhedron.getViewName(index);
	}

	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#getTransformation()
	 */
	@Override
	public Matrix4d getTransformation() {	
	    return rotationAxisAligner.getTransformation();
	}
	
	public void setDefaultColoring(String colorScript) {
		this.defaultColoring = colorScript;
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#drawPolyhedron()
	 */
	@Override
	public String drawPolyhedron() {
		StringBuilder s = new StringBuilder();

		Point3d[] vertices = getPolyhedronVertices();
		
		int index = 0;
		double width = getMaxExtension()*0.015;

		for (int[] lineLoop: polyhedron.getLineLoops()) {
			s.append("draw polyhedron");
			s.append(name);
			s.append(index++);
			s.append(" line");
			for (int i: lineLoop) {
				s.append(getJmolPoint(vertices[i]));
			}
			s.append("width ");
		    s.append(fDot2(width));
			s.append(" color");
			Color4f c = getPolyhedronColor();
			s.append(getJmolColor(c));
			s.append(" off;");
		}

		return s.toString();
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#hidePolyhedron()
	 */
	@Override
	public String hidePolyhedron() {
		return "draw polyhedron* off;";
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#showPolyhedron()
	 */
	@Override
	public String showPolyhedron() {
		return "draw polyhedron* on;";
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#drawAxes()
	 */
	@Override
	public String drawAxes() {
		if (rotationGroup.getPointGroup().equals("C1")) {
			return drawInertiaAxes();
		} else {
			return drawSymmetryAxes();
		}
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#hideAxes()
	 */
	@Override
	public String hideAxes() {
		return "draw axes* off;";
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#showAxes()
	 */
	@Override
	public String showAxes() {
		return "draw axes* on;";
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#playOrientations()
	 */
	@Override
	public String playOrientations() {
		StringBuilder s = new StringBuilder();
		
		// draw point group
				
		if ( rotationGroup.getPointGroup().equals("C1")) {
			s.append(drawFooter("Asymmetric", "white"));
		} else {
			s.append(drawFooter("Point group " + rotationGroup.getPointGroup(), "white"));
		}
		
		// draw polygon
		s.append(drawPolyhedron()); // draw invisibly
		s.append(showPolyhedron());
			
		// draw axes
		s.append(drawAxes());
		s.append(showAxes());
		
		// loop over all orientations with 4 sec. delay
		for (int i = 0; i < getOrientationCount(); i++) {
			s.append(deleteHeader());
			s.append(getOrientationWithZoom(i));
			s.append(drawHeader(polyhedron.getViewName(i), "white"));
			s.append("delay 4;");
		}
		
		// go back to first orientation
		s.append(deleteHeader());
		s.append(getOrientationWithZoom(0));
		s.append(drawHeader(polyhedron.getViewName(0), "white"));
		
		return s.toString();
	}	
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#colorBySubunit()
	 */
	@Override
	public String colorBySubunit() {
	    Subunits subunits = rotationAxisAligner.getSubunits();
	    List<Integer> modelNumbers = subunits.getModelNumbers();
	    List<String> chainIds = subunits.getChainIds();
	    List<List<Integer>> orbits = rotationAxisAligner.getOrbits();
		int fold = rotationGroup.getRotation(0).getFold();

		Color[] col = null;
		Color4f[] colors = null;
		if (fold > 1) {
	        col = ColorBrewer.Spectral.getColorPalette(2*fold);
	        colors = ColorConverter.convertColor4f(col);
		} else {
			col = ColorBrewer.Spectral.getColorPalette(orbits.size());
	        colors = ColorConverter.convertColor4f(col);
		}
		int half = colors.length/2;
		for (int i = 0; i < half; i++) {
			if (i % 2 != 0) {
			   Color4f temp = colors[i];
			   colors[i] = colors[half+i];
			   colors[half+i] = temp;
			}
		}
	    Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();
	    
		for (int i = 0; i < orbits.size(); i++) {
			for (int j = 0; j < fold; j++) {
				// assign alternating color sets to adjacent orbits
				int colorIndex = i;
				if (fold > 1) {
					if (i % 2 == 0) {
						colorIndex = j;
					} else {
						colorIndex = fold + j;
					}
				}
				int subunit = orbits.get(i).get(j);
				Color4f c = colors[colorIndex];
				List<String> ids = colorMap.get(c);
				if (ids == null) {
					ids = new ArrayList<String>();
					colorMap.put(c,  ids);
				}
				String id = getChainSpecification(modelNumbers, chainIds, subunit);
				ids.add(id);
			}
		}
		return defaultColoring + getJmolColorScript(colorMap) + getJmolLigandScript();
	}

	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#colorBySequenceCluster()
	 */
	@Override
	public String colorBySequenceCluster() {
	    Subunits subunits = rotationAxisAligner.getSubunits();
	    int n = subunits.getSubunitCount();
	    List<Integer> modelNumbers = subunits.getModelNumbers();
	    List<String> chainIds = subunits.getChainIds();
	    List<Integer> seqClusterIds = subunits.getSequenceClusterIds();
	    int clusters = Collections.max(seqClusterIds) + 1;
	    Color[] col = ColorBrewer.BrBG.getColorPalette(clusters);
	    Color4f[] colors = ColorConverter.convertColor4f(col);
		Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();
		
		for (int i = 0; i < n; i++) {
			Color4f c = colors[seqClusterIds.get(i)];
			List<String> ids = colorMap.get(c);
			if (ids == null) {
				ids = new ArrayList<String>();
				colorMap.put(c,  ids);
			}
			String id = getChainSpecification(modelNumbers, chainIds, i);
			ids.add(id);

		}
		
		return defaultColoring + getJmolColorScript(colorMap) + getJmolLigandScript();
	}
	
	/* (non-Javadoc)
	 * @see org.biojava.nbio.structure.quaternary.jmolScript.JMolSymmetryScriptInterface#colorBySymmetry()
	 */
	@Override
	public String colorBySymmetry() {
		// TODO needs some refactoring
		String pointGroup = rotationGroup.getPointGroup();
		Subunits subunits = rotationAxisAligner.getSubunits();
		List<Integer> modelNumbers = subunits.getModelNumbers();
		List<String> chainIds = subunits.getChainIds();
		List<List<Integer>> orbits = rotationAxisAligner.getOrbits();

		int n = subunits.getSubunitCount();
		int fold = rotationGroup.getRotation(0).getFold();
		
		Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();

		// Simple Cn symmetry
		if (pointGroup.startsWith("C") && n == fold) {
			colorMap = getCnColorMap();
			// complex cases
		} else if ((pointGroup.startsWith("D") && orbits.size() > 2) || 
				pointGroup.equals("T")|| pointGroup.equals("O") || pointGroup.equals("I")) {
			int nColor = 0;
			if (orbits.size() % 2 == 0) {
				nColor = orbits.size()/2;
			} else {
				nColor = (orbits.size() + 1)/2;
			}
			Color4f[] colors = getSymmetryColors(nColor); 

			for (int i = 0; i < orbits.size(); i++) {
				int colorIndex = i;
				
				// reverse colors once the center of the structure has been reached
				if (i >= nColor) {
					colorIndex = orbits.size() - 1 - i;
				}
				Color4f c = colors[colorIndex];
				List<String> ids = colorMap.get(c);
				if (ids == null) {
					ids = new ArrayList<String>();
					colorMap.put(c,  ids);
				}
				for (int subunit: orbits.get(i)) {
					String id = getChainSpecification(modelNumbers, chainIds, subunit);
					ids.add(id);
				}
			}

			// Simple Dn symmetry
		} else {
			Color4f[] colors = getSymmetryColors(orbits.size());
			
			for (int i = 0; i < orbits.size(); i++) {
				Color4f c = new Color4f(colors[i]);
				List<String> ids = colorMap.get(c);
				if (ids == null) {
					ids = new ArrayList<String>();
					colorMap.put(c,  ids);
				}
				List<Integer> orbit = orbits.get(i);
				for (int j = 0; j < orbit.size(); j++) {
					String id = getChainSpecification(modelNumbers, chainIds, orbit.get(j));
					ids.add(id);
				}
			}
		}
		
		return defaultColoring + getJmolColorScript(colorMap) + getJmolLigandScript();
	}
	
	// --- protected methods ---
	/**
	 * Returns the maximum extension (length) of structure
	 * @return
	 */
	protected double getMaxExtension() {
		Vector3d dimension = rotationAxisAligner.getDimension();
		double maxExtension = Math.max(dimension.x, dimension.y);
		maxExtension = Math.max(maxExtension, dimension.z);
		return maxExtension;
	}
	
	/**
	 * Returns the mean extension (length) of structure
	 * @return
	 */
	protected double getMeanExtension() {
		Vector3d dimension = rotationAxisAligner.getDimension();
		return (dimension.x+dimension.y+dimension.z)/3;
	}
	
	/**
	 * @return the axisTransformation
	 */
	protected RotationAxisAligner getAxisTransformation() {
		return rotationAxisAligner;
	}

	/**
	 * @param axisTransformation the axisTransformation to set
	 */
	protected void setAxisTransformation(RotationAxisAligner axisTransformation) {
		this.rotationAxisAligner = axisTransformation;
	}

	/**
	 * @return the rotationGroup
	 */
	protected RotationGroup getRotationGroup() {
		return rotationGroup;
	}

	/**
	 * @param rotationGroup the rotationGroup to set
	 */
	protected void setRotationGroup(RotationGroup rotationGroup) {
		this.rotationGroup = rotationGroup;
	}

	/**
	 * @return the polyhedron
	 */
	protected Polyhedron getPolyhedron() {
		return polyhedron;
	}

	/**
	 * @param polyhedron the polyhedron to set
	 */
	protected void setPolyhedron(Polyhedron polyhedron) {
		this.polyhedron = polyhedron;
	}
	
//  --- private methods ---
	
	private String getChainSpecification(List<Integer> modelNumbers, List<String> chainIds, int subunit) {
		if (onTheFly) {
			if (Collections.max(modelNumbers) > 1) {
				return chainIds.get(subunit) + "&symop=" + (modelNumbers.get(subunit)+1);
			} else {
				// if there is only a single symop, Jmol does not accept the symop syntax
				return chainIds.get(subunit);
			}
		} else {
			return chainIds.get(subunit) + "/" + (modelNumbers.get(subunit)+1);
		}
	}
	
	private Map<Color4f, List<String>> getCnColorMap() {
		Subunits subunits = rotationAxisAligner.getSubunits();
		List<Integer> modelNumbers = subunits.getModelNumbers();
		List<String> chainIds = subunits.getChainIds();
		List<List<Integer>> orbits = rotationAxisAligner.getOrbits();

		int fold = rotationGroup.getRotation(0).getFold();

		Map<Color4f, List<String>> colorMap = new HashMap<Color4f, List<String>>();
		Color4f[] colors = getSymmetryColors(fold);

		for (List<Integer> orbit: orbits) {
			for (int i = 0; i < fold; i++) {
				int subunit = orbit.get(i);
				String id = null;
				id = getChainSpecification(modelNumbers, chainIds, subunit);
				Color4f c = colors[i];
				List<String> ids = colorMap.get(c);
				if (ids == null) {
					ids = new ArrayList<String>();
					colorMap.put(c, ids);
				}
				ids.add(id);
			}
		}

		return colorMap;
	}
	
	private Point3d[] getPolyhedronVertices() {
		Point3d[] vertices = polyhedron.getVertices();
		Matrix4d reverseTransformation = rotationAxisAligner.getGeometicCenterTransformation();
		for (int i = 0; i < vertices.length; i++) {
			reverseTransformation.transform(vertices[i]);
		}
		return vertices;
	}
	
	/**
	 * Return a color that is complementary to the symmetry color
	 * @return
	 */
	private Color4f getPolyhedronColor() {
		Color4f[] colors = getSymmetryColors(5);
		Color4f strongestColor = colors[4];
		Color4f complement = new Color4f(Color.WHITE);
		complement.sub(strongestColor);
		return complement;
	}
	
	/**
	 * Returns a unique color palette based on point group
	 * @param nColors
	 * @return
	 */
	private Color4f[] getSymmetryColors(int nColors) {
		String pointGroup = rotationGroup.getPointGroup();
		Color[] col = null;
		Color4f[] colors = null;
		if (pointGroup.equals("C1")) {
			col = ColorBrewer.Greys.getColorPalette(nColors);
	        colors = ColorConverter.convertColor4f(col);
		} else if (pointGroup.startsWith("C")) {
			col = ColorBrewer.YlGnBu.getColorPalette(nColors);	
	        colors = ColorConverter.convertColor4f(col);
		} else if (pointGroup.startsWith("D")) {
			col = ColorBrewer.YlOrRd.getColorPalette(nColors);
	        colors = ColorConverter.convertColor4f(col);
		} else if (pointGroup.equals("T")) {
			col = ColorBrewer.Greens.getColorPalette(nColors);
	        colors = ColorConverter.convertColor4f(col);
		} else if (pointGroup.equals("O")) {
			col = ColorBrewer.Blues.getColorPalette(nColors);
	        colors = ColorConverter.convertColor4f(col);
		} else if (pointGroup.equals("I")) {
			col = ColorBrewer.BuPu.getColorPalette(nColors);
	        colors = ColorConverter.convertColor4f(col);
		} else {
			col = ColorBrewer.Greys.getColorPalette(nColors);
	        colors = ColorConverter.convertColor4f(col);
		}
		System.arraycopy(colors, 0, colors, 0, nColors);
		return colors;	
	}
	
	private String drawInertiaAxes() {
		StringBuilder s = new StringBuilder();
		Point3d centroid = rotationAxisAligner.getGeometricCenter();
		Vector3d[] axes = rotationAxisAligner.getPrincipalAxesOfInertia();

		for (int i = 0; i < axes.length; i++) {
			s.append("draw axesInertia");
			s.append(name);
			s.append(i);
			s.append(" ");
			s.append("line");
			Point3d v1 = new Point3d(axes[i]);
			if (i == 0) {
				v1.scale(AXIS_SCALE_FACTOR*rotationAxisAligner.getDimension().y);
			} else if (i == 1) {
				v1.scale(AXIS_SCALE_FACTOR*rotationAxisAligner.getDimension().x);
			} else if (i == 2) {
				v1.scale(AXIS_SCALE_FACTOR*rotationAxisAligner.getDimension().z);
			}
			Point3d v2 = new Point3d(v1);
			v2.negate();
			v1.add(centroid);
			v2.add(centroid);
			s.append(getJmolPoint(v1));
			s.append(getJmolPoint(v2));
			s.append("width 0.5 ");
			s.append(" color white");
			s.append(" off;");
		}
        return s.toString();
	};
	
	private String drawSymmetryAxes() {
		StringBuilder s = new StringBuilder();

		int n = rotationGroup.getOrder();
		if (n == 0) {
			return s.toString();
		}

		float diameter = 0.5f;
		double radius = 0;
		String color = "";

		List<Rotation> axes = getUniqueAxes();

		int i = 0;
		for (Rotation r: axes) {
			if (rotationGroup.getPointGroup().startsWith("C") || (rotationGroup.getPointGroup().startsWith("D") && r.getDirection() == 0)) {
				radius =  rotationAxisAligner.getDimension().z; // principal axis uses z-dimension
				color = N_FOLD_AXIS_COLOR;
			} else {
				radius = polyhedron.getCirumscribedRadius();
			
				if (r.getFold() == 2) {
					color = TWO_FOLD_AXIS_COLOR;
				} else if (r.getFold() == 3) {
					color = THREE_FOLD_AXIS_COLOR;
				} else {
					color = N_FOLD_AXIS_COLOR;
				}
			}
		

			Point3d center = rotationAxisAligner.getGeometricCenter();
			AxisAngle4d axisAngle = r.getAxisAngle();
			Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
			Vector3d refAxis = rotationAxisAligner.getRotationReferenceAxis();
			
			s.append(getSymmetryAxis(i, i+axes.size(), rotationGroup.getPointGroup(), r.getFold(), refAxis, radius, diameter, color, center, axis));
	        i++;
		}

		return s.toString();
	}

	private Vector3d getAligmentVector(Point3d point, Vector3d axis) {		
		// for system with a single Cn axis
		if (rotationGroup.getPointGroup().startsWith("C") || rotationGroup.getPointGroup().startsWith("D")) {
			// if axis is orthogonal to principal axis, use principal axis as reference axis
			if (axis.dot(rotationAxisAligner.getPrincipalRotationAxis()) < 0.1) {
				return rotationAxisAligner.getPrincipalRotationAxis();
			} else {
				return rotationAxisAligner.getRotationReferenceAxis();
			}
		}

		// for T, O, and I point groups find reference axis with respect to
		// nearest polyhedron vertex
		Vector3d ref = new Vector3d();
		double dSqThreshold = 25;
		double dSqMin = Double.MAX_VALUE;
		Point3d refPoint = null;
		// find nearest point on polyhedron as reference point,
		// but do not choose a point on the same axis (therefore, we 
		// apply a distance threshold squared 5A*5A = 25A^2
		for (Point3d v: getPolyhedronVertices()) {
			double dSq = point.distanceSquared(v);
			if (dSq > dSqThreshold) {
				if (dSq < dSqMin) {
					dSqMin = dSq;
					refPoint = v;
				}
			}
		}


		ref.sub(point, refPoint);

		// this ref vector is usually not orthogonal to the 
		// axis. The following double-cross product makes it
		// orthogonal.
		ref.cross(axis, ref);
		ref.cross(axis, ref); // note, done twice on purpose
		ref.normalize();
		return ref;
	}
	
	private String getSymmetryAxis(int i, int j, String pointGroup, int n, Vector3d referenceAxis, double radius, float diameter, String color, Point3d center, Vector3d axis) {
		boolean drawPolygon = true;
		
		Point3d p1 = new Point3d(axis);
		p1.scaleAdd(-AXIS_SCALE_FACTOR * radius, center);

		Point3d p2 = new Point3d(axis);
		p2.scaleAdd(AXIS_SCALE_FACTOR * radius, center);

		StringBuilder s = new StringBuilder();
		s.append("draw");
		s.append(" axesSymmetry");
		s.append(name);
		s.append(i);
		s.append(" cylinder");
		s.append(getJmolPoint(p1));
		s.append(getJmolPoint(p2));
		s.append("diameter ");
		s.append(diameter);
		s.append(" color ");
		s.append(color);
		s.append(" off;");

		// calc. point to center symmetry symbols. They are offset by 0.01
		// to avoid overlap with the polyhedron
		p1 = new Point3d(axis);
		p1.scaleAdd(-1.01*radius, center);

		p2 = new Point3d(axis);
		p2.scaleAdd(1.01*radius, center);

		if (drawPolygon == true) {
			double polygonRadius = getMeanExtension() * 0.06;
			if (n == 2) {
				referenceAxis = getAligmentVector(p1, axis);
				s.append(getC2PolygonJmol(i, p1, referenceAxis, axis, color, polygonRadius, name));
				referenceAxis = getAligmentVector(p2, axis);
				s.append(getC2PolygonJmol(j, p2,  referenceAxis, axis, color, polygonRadius, name));
			} else if (n > 2) {
				referenceAxis = getAligmentVector(p1, axis);
				s.append(getPolygonJmol(i, p1, referenceAxis, axis, n, color, polygonRadius, name));
				referenceAxis = getAligmentVector(p2, axis);
				s.append(getPolygonJmol(j, p2, referenceAxis, axis, n, color, polygonRadius, name));
			}
		}

		return s.toString();
	}
	
	private static String getPolygonJmol(int index, Point3d center, Vector3d referenceAxis, Vector3d axis, int n, String color, double radius, String name) {
		StringBuilder s = new StringBuilder();
		s.append("draw axesSymbol");
		s.append(name);
		s.append(index);
		s.append(" ");
		s.append("polygon");
		s.append(" ");
		s.append(n+1); 
		s.append(getJmolPoint(center));

		Vector3d[] vertexes = getPolygonVertices(axis, referenceAxis, center, n, radius);
		// create vertex list
		for (Vector3d v: vertexes) {
			s.append(getJmolPoint(v));
		}

		// create face list
		s.append(n);
		for (int i = 1; i <= n; i++) {
			s.append("[");
			s.append(0);
			s.append(" ");
			s.append(i);
			s.append(" ");
			if (i < n) {
				s.append(i+1);
			} else {
				s.append(1);
			}
			s.append(" ");
			s.append(7);
			s.append("]");
		}

		if (n == 2) {
	      	s.append("mesh off");
		}
		s.append(" color ");
		s.append(color);
		s.append(" off;");

		return s.toString();
	}
	
	private static Vector3d[] getPolygonVertices(Vector3d axis, Vector3d referenceAxis, Point3d center, int n, double radius) {
		Vector3d ref = new Vector3d(referenceAxis);
		ref.scale(radius);		

		AxisAngle4d axisAngle = new AxisAngle4d(axis, 0);
		Vector3d[] vectors = new Vector3d[n];
		Matrix4d m = new Matrix4d();

		for (int i = 0; i < n; i++) {
			axisAngle.angle = i * 2 * Math.PI/n;
			vectors[i] = new Vector3d(ref);		
			m.set(axisAngle);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m.setElement(3, 3, 1.0);
			m.transform(vectors[i]);
			vectors[i].add(center);
		}
		return vectors;
	}
	
	private static String getC2PolygonJmol(int index, Point3d center, Vector3d referenceAxis, Vector3d axis, String color, double radius, String name) {
		StringBuilder s = new StringBuilder();
		int n = 10;
		s.append("draw axesSymbol");
		s.append(name);
		s.append(index);
		s.append(" ");
		s.append("polygon");
		s.append(" ");
		s.append(n-1); 
		s.append(getJmolPoint(center));

		Vector3d[] vertexes = getC2PolygonVertices(axis, referenceAxis, center, n, radius);
		// create vertex list
		for (Vector3d v: vertexes) {
			s.append(getJmolPoint(v));
		}

		// create face list
		s.append(n-2);

		for (int i = 1; i < n-1; i++) {
			s.append("[");
			s.append(0);
			s.append(" ");
			s.append(i);
			s.append(" ");
			if (i < n-2) {
				s.append(i+1);
			} else {
				s.append(1);
			}
			s.append(" ");
			s.append(7);
			s.append("]");
		}

		s.append("color ");
		s.append(color);
		s.append(" off;");

		return s.toString();
	}
	private static Vector3d[] getC2PolygonVertices(Vector3d axis, Vector3d referenceAxis, Point3d center, int n, double radius) {
		Vector3d ref = new Vector3d(referenceAxis);
		ref.scale(4*radius);		

		AxisAngle4d axisAngle = new AxisAngle4d(axis, 0);
		int k = n / 2;
		// draw arc 1/6 of a full circle
		int f = 6;
		Vector3d[] vectors = new Vector3d[n-2];	
		Matrix4d m = new Matrix4d();
		
		// first point of arc
		axisAngle.angle = (k+0.5) * 2 * Math.PI/(f*k);
		Vector3d begin = new Vector3d(ref);		
		m.set(axisAngle);
		// make sure matrix element m33 is 1.0. It's 0 on Linux.
		m.setElement(3, 3, 1.0);
		m.transform(begin);
		
		// last point of arc
		axisAngle.angle = (2*k-1+0.5) * 2 * Math.PI/(f*k);
		Vector3d end = new Vector3d(ref);		
		m.set(axisAngle);
		// make sure matrix element m33 is 1.0. It's 0 on Linux.
		m.setElement(3, 3, 1.0);
		m.transform(end);
		
		// center of arc
		Vector3d arcCenter = new Vector3d();
		arcCenter.interpolate(begin, end, 0.5);
		arcCenter.negate();
		
		// add translation component
		Vector3d trans =  new Vector3d();
		trans.sub(center, arcCenter);

		// draw arc (1/6 of a full circle)
		for (int i = 0; i < k; i++) {
			axisAngle.angle = (k + i + 0.5) * 2 * Math.PI/(f*k);
			vectors[i] = new Vector3d(ref);		
			m.set(axisAngle);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m.setElement(3, 3, 1.0);
			m.transform(vectors[i]);
			vectors[i].add(arcCenter);
			vectors[i].add(center);
		}
		// in reverse order, draw reflected half of arc (1/6 of full circle)
		// don't draw first and last element, since the are already part of the previous arc
		for (int i = k; i < 2*k-2; i++) {
			axisAngle.angle = (f/2*k + i + 1.5) * 2 * Math.PI/(f*k);
			vectors[i] = new Vector3d(ref);		
			m.set(axisAngle);
			// make sure matrix element m33 is 1.0. It's 0 on Linux.
			m.setElement(3, 3, 1.0);
			m.transform(vectors[i]);
			vectors[i].sub(arcCenter);
			vectors[i].add(center);
		}
		return vectors;
	}
	

    private List<Rotation> getUniqueAxes() {
    	List<Rotation> uniqueRotations = new ArrayList<Rotation>();
    	
    	for (int i = 0, n = rotationGroup.getOrder(); i < n; i++) {
			Rotation rotationI = rotationGroup.getRotation(i);
			AxisAngle4d axisAngleI = rotationI.getAxisAngle();
			Vector3d axisI = new Vector3d(axisAngleI.x, axisAngleI.y, axisAngleI.z);
			
			boolean redundant = false;
			for (Rotation r: uniqueRotations) {
				AxisAngle4d axisAngleJ = r.getAxisAngle();
			    Vector3d axisJ = new Vector3d(axisAngleJ.x, axisAngleJ.y, axisAngleJ.z);
			    if (Math.abs(axisI.dot(axisJ)) > 0.99) {
			    	redundant = true;
			    	break;
			    }
			}
			
			if (! redundant) {
				uniqueRotations.add(rotationI);
			}
    	}
        return uniqueRotations;
    }
	
	private String drawHeader(String text, String color) {
		StringBuilder s = new StringBuilder();
		s.append("set echo top center;");
		s.append("color echo ");
		s.append(color);
		s.append(";");
		s.append("font echo 24 sanserif;");
		s.append("echo ");
		s.append(text);
		s.append(";");
		return s.toString();
	}
	
	private String deleteHeader() {
		return "set echo top center;echo ;";
	}
	
	private String drawFooter(String text, String color) {
		StringBuilder s = new StringBuilder();
		s.append("set echo bottom center;");
		s.append("color echo ");
		s.append(color);
		s.append(";");
		s.append("font echo 24 sanserif;");
		s.append("echo "+ text);
		//s.append("echo Point group ");
		//s.append(rotationGroup.getPointGroup());
		s.append(";");
		return s.toString();
	}
	
	private String setCentroid() {
		// calculate center of rotation
	//	Point3d centroid = axisTransformation.getGeometricCenter();
		Point3d centroid = rotationAxisAligner.getCentroid();
			
		// set centroid
		StringBuilder s = new StringBuilder();
		s.append("center");
		s.append(getJmolPoint(centroid));
		s.append(";");
		return s.toString();
	}
}
