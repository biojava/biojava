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
package org.biojava.nbio.structure.xtal;

import org.biojava.nbio.structure.jama.EigenvalueDecomposition;
import org.biojava.nbio.structure.jama.Matrix;
import org.biojava.nbio.structure.xtal.io.TransfAlgebraicAdapter;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * A crystallographic space group. We store the standard numeric identifier,
 * the international short symbol and the transformations corresponding to
 * each space group (as Matrix4ds and in algebraic notation).
 * The information for all (protein crystallography) space groups can be
 * parsed from the XML file in the resource directory.
 *
 * See: http://en.wikipedia.org/wiki/Space_group
 *
 * @author duarte_j
 * @see SymoplibParser
 */
@XmlRootElement(name = "SpaceGroup", namespace ="http://www.biojava.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class SpaceGroup implements Serializable {

	private static final long serialVersionUID = 1L;
	private static final Logger logger = LoggerFactory.getLogger(SpaceGroup.class);


	private static final Pattern splitPat1 = Pattern.compile("((?:[+-]?[XYZ])+)([+-][0-9/.]+)");
	private static final Pattern splitPat2 = Pattern.compile("([+-]?[0-9/.]+)((?:[+-][XYZ])+)");
	private static final Pattern coordPat = Pattern.compile("(?:([+-])?([XYZ]))+?"); // the last +? is for ungreedy matching
	private static final Pattern transCoefPat = Pattern.compile("([-+]?[0-9.]+)(?:/([0-9.]+))?");

	private static final Pattern nonEnantPat = Pattern.compile("[-abcmnd]");

	protected static final double DELTA=0.0000001;

	private  int id;
	private  int multiplicity;
	private  int primitiveMultiplicity;
	private  String shortSymbol;
	private  String altShortSymbol;
	private  List<Matrix4d> transformations;
	private  List<String> transfAlgebraic;
	private  Vector3d[] cellTranslations; // in space groups I, C, F or H there are pure cell translations corresponding to recenterings

	private AxisAngle4d[] axisAngles;

	private int[] axisTypes; // indices of array are transformIds

	private BravaisLattice bravLattice;

	@SuppressWarnings("unused")
	private SpaceGroup(){
		// required by JAXB

	}

	public SpaceGroup(int id, int multiplicity, int primitiveMultiplicity, String shortSymbol, String altShortSymbol, BravaisLattice bravLattice) {
		this.id = id;
		this.multiplicity = multiplicity;
		this.primitiveMultiplicity = primitiveMultiplicity;
		this.shortSymbol = shortSymbol;
		this.altShortSymbol = altShortSymbol;
		transformations = new ArrayList<Matrix4d>(multiplicity);
		transfAlgebraic = new ArrayList<String>(multiplicity);
		cellTranslations = new Vector3d[multiplicity/primitiveMultiplicity];
		this.bravLattice = bravLattice;
	}

	/**
	 * Get the space group for the given international short name, using
	 * the PDB format, e.g. 'P 21 21 21' or 'C 1 c 1'
	 * @param shortName
	 * @return the SpaceGroup or null if the shortName is not valid
	 * @see SymoplibParser#getSpaceGroup(String)
	 */
	public static SpaceGroup parseSpaceGroup(String shortName) {
		return SymoplibParser.getSpaceGroup(shortName);
	}

	public void addTransformation(String transfAlgebraic) {
		this.transfAlgebraic.add(transfAlgebraic);
		this.transformations.add(getMatrixFromAlgebraic(transfAlgebraic));
	}

	protected void initializeCellTranslations() {
		if ( cellTranslations != null && cellTranslations.length >0) {
			// we already initialized this
			return;
		}
		cellTranslations = new Vector3d[multiplicity/primitiveMultiplicity];
		cellTranslations[0] = new Vector3d(0,0,0);

		if ( transformations == null){
			logger.warn("transformations == null" + this.toXML());
		}

		if (multiplicity==primitiveMultiplicity) {
			return;
		}
		int fold = multiplicity/primitiveMultiplicity;



		for (int n=1;n<fold;n++) {
			if ( transformations.size() < (n* primitiveMultiplicity)){
				logger.warn("WARNING number of transformations < " +(n*primitiveMultiplicity));
				logger.warn(this.toXML());
			}
			Matrix4d t = transformations.get(n*primitiveMultiplicity);
			cellTranslations[n] = new Vector3d(t.m03,t.m13,t.m23);
		}
	}

	public int getMultiplicity() {
		return multiplicity;
	}

	public int getPrimitiveMultiplicity() {
		return primitiveMultiplicity;
	}

	public Vector3d[] getCellTranslations() {
		return cellTranslations;
	}

	public Vector3d getCellTranslation(int i) {
		return cellTranslations[i];
	}

	public static Matrix4d getMatrixFromAlgebraic(String transfAlgebraic) {
		String[] parts = transfAlgebraic.toUpperCase().split(",");
		double[] xCoef = convertAlgebraicStrToCoefficients(parts[0].trim());
		double[] yCoef = convertAlgebraicStrToCoefficients(parts[1].trim());
		double[] zCoef = convertAlgebraicStrToCoefficients(parts[2].trim());

		Matrix4d mat = new Matrix4d();
		mat.setIdentity();
		mat.setRotation(new Matrix3d(xCoef[0],xCoef[1],xCoef[2],yCoef[0],yCoef[1],yCoef[2],zCoef[0],zCoef[1],zCoef[2]));
		mat.setTranslation(new Vector3d(xCoef[3],yCoef[3],zCoef[3]));
		return mat;
		//return new Matrix4d(xCoef[0],xCoef[1],xCoef[2],xCoef[3],
		//					yCoef[0],yCoef[1],yCoef[2],yCoef[3],
		//					zCoef[0],zCoef[1],zCoef[2],zCoef[3],
		//					0,0,0,1);
	}

	private static double[] convertAlgebraicStrToCoefficients(String algString) {
		String letters = null;
		String noLetters = null;
		Matcher m = splitPat1.matcher(algString);
		if (m.matches()) {
			letters = m.group(1);
			noLetters = m.group(2);
		} else {
			m = splitPat2.matcher(algString);
			if (m.matches()) {
				letters = m.group(2);
				noLetters = m.group(1);
			} else {
				letters = algString;
			}
		}
		double[] coefficients = new double[4];
		m = coordPat.matcher(letters);
		while(m.find()){
			String sign = "";
			if (m.group(1)!=null) {
				sign = m.group(1);
			}
			double s = 1.0;
			if (sign.equals("-")){
				s = -1.0;
			}
			String coord = m.group(2);
			if (coord.equals("X")) {
				coefficients[0] = s;
			} else if (coord.equals("Y")) {
				coefficients[1] = s;
			} else if (coord.equals("Z")) {
				coefficients[2] = s;
			}
		}
		if (noLetters!=null) {
			m = transCoefPat.matcher(noLetters);
			if (m.matches()) {
				double num = Double.parseDouble(m.group(1));
				double den = 1;
				if (m.group(2)!=null) {
					den = Double.parseDouble(m.group(2));
				}
				coefficients[3] = num/den;
			}
		} else {
			coefficients[3]=0;
		}
		return coefficients;
	}

	/**
	 * Gets the standard numeric identifier for the space group.
	 * See for example http://en.wikipedia.org/wiki/Space_group
	 * or the IUCr crystallographic tables
	 * @return
	 */
	public int getId() {
		return id;
	}

	/**
	 * Gets the international short name (as used in PDB),
	 * e.g. "P 21 21 21" or "C 1 c 1"
	 * @return
	 */
	public String getShortSymbol() {
		return shortSymbol;
	}

	/**
	 * Gets the alternative international short name (as sometimes used in PDB),
	 * e.g. "I 1 2 1" instead of "I 2"
	 * @return
	 */
	public String getAltShortSymbol() {
		return altShortSymbol;
	}

	/**
	 * Gets all transformations except for the identity in crystal axes basis.
	 * @return
	 */
	public List<Matrix4d> getTransformations() {
		List<Matrix4d> transfs = new ArrayList<Matrix4d>();
		for (int i=1;i<this.transformations.size();i++){
			transfs.add(transformations.get(i));
		}
		return transfs;
	}

	private void calcRotAxesAndAngles() {

		axisAngles = new AxisAngle4d[multiplicity];

		// identity operator (transformId==0)
		axisAngles[0] = new AxisAngle4d(new Vector3d(0,0,0), 0.0);

		for (int i=1;i<this.transformations.size();i++){
			Matrix3d r = new Matrix3d(transformations.get(i).m00,transformations.get(i).m01,transformations.get(i).m02,
					transformations.get(i).m10,transformations.get(i).m11,transformations.get(i).m12,
					transformations.get(i).m20,transformations.get(i).m21,transformations.get(i).m22);

			axisAngles[i] = getRotAxisAndAngle(r);

		}
	}

	/**
	 * Calculates the axis fold type (1, 2, 3, 4, 5, 6 for rotations or -1, -2, -3, -4, -6 improper rotations)
	 * from the trace of the rotation matrix, see for instance
	 * http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf
	 */
	private void calcAxisFoldTypes() {
		axisTypes = new int[multiplicity];

		for (int i=0;i<this.transformations.size();i++){

			axisTypes[i] = getRotAxisType(transformations.get(i));

		}
	}

	public AxisAngle4d getRotAxisAngle(int transformId) {
		if (this.axisAngles == null) calcRotAxesAndAngles();
		return this.axisAngles[transformId];
	}

	/**
	 * Returns true if both given transform ids belong to the same crystallographic axis (a, b or c)
	 * For two non-rotation transformations (i.e. identity operators) it returns true
	 * @param tId1
	 * @param tId2
	 * @return
	 */
	public boolean areInSameAxis(int tId1, int tId2) {
		if (tId1==tId2) return true;

		if (axisAngles== null) calcRotAxesAndAngles();

		if (getAxisFoldType(tId1)==1 && getAxisFoldType(tId2)==1) return true;

		// we can't deal yet with improper rotations: we return false whenever either of them is improper
		if (getAxisFoldType(tId1)<0 || getAxisFoldType(tId2)<0) return false;

		Vector3d axis1 = new Vector3d(axisAngles[tId1].x, axisAngles[tId1].y, axisAngles[tId1].z);
		Vector3d axis2 = new Vector3d(axisAngles[tId2].x, axisAngles[tId2].y, axisAngles[tId2].z);

		// TODO revise: we might need to consider that the 2 are in same direction but opposite senses
		// the method is not used at the moment anyway
		if (deltaComp(axis1.angle(axis2), 0.0, DELTA)) return true;

		return false;
	}

	/**
	 * Given a transformId returns the type of axis of rotation: 1 (no rotation), 2, 3, 4 or 6 -fold
	 * and for improper rotations: -1, -2, -3, -4 and -6
	 *
	 * @param transformId
	 * @return
	 */
	public int getAxisFoldType(int transformId) {
		if (axisTypes== null) calcAxisFoldTypes();
		return axisTypes[transformId];
	}

	/**
	 * Gets a transformation by index expressed in crystal axes basis.
	 * Index 0 corresponds always to the identity transformation.
	 * Beware the returned Matrix4d is not a copy but it stays linked
	 * to the one stored in this SpaceGroup object
	 * @param i
	 * @return
	 */
	public Matrix4d getTransformation(int i) {
		return transformations.get(i);
	}

	/**
	 * Gets a transformation algebraic string given its index.
	 * Index 0 corresponds always to the identity transformation.
	 * @param i
	 * @return
	 */
	public String getTransfAlgebraic(int i) {
		return transfAlgebraic.get(i);
	}

	@Override
	public int hashCode() {
	    final int prime = 31;
	    int result = 1;
	    result = prime * result + id;
	    return result;
	}

	@Override
	public boolean equals(Object o) {
		if (! (o instanceof SpaceGroup)) {
			return false;
		}
		SpaceGroup other = (SpaceGroup) o;
		if (other.getId()==this.getId()) {
			return true;
		}
		return false;
	}

	/**
	 * Gets the number of symmetry operators corresponding to this SpaceGroup (counting
	 * the identity operator)
	 * @return
	 */
	public int getNumOperators() {
		return this.transformations.size();
	}

	public static String getAlgebraicFromMatrix(Matrix4d m) {
		String x = formatAlg(m.m00,m.m01,m.m02,m.m03);
		String y = formatAlg(m.m10,m.m11,m.m12,m.m13);
		String z = formatAlg(m.m20,m.m21,m.m22,m.m23);
		String alg = x+","+y+","+z;
		return alg;
	}

	private static String formatAlg(double xcoef, double ycoef, double zcoef, double trans) {
		boolean[] leading = {false,false,false};
		if (xcoef!=0) {
			leading[0] = true;
		} else if (ycoef!=0) {
			leading[1] = true;
		} else if (zcoef!=0) {
			leading[2] = true;
		}
		String x = deltaComp(xcoef,0,DELTA)?"":formatCoef(xcoef,leading[0])+"X";
		String y = deltaComp(ycoef,0,DELTA)?"":formatCoef(ycoef,leading[1])+"Y";
		String z = deltaComp(zcoef,0,DELTA)?"":formatCoef(zcoef, leading[2])+"Z";
		String t = deltaComp(trans,0,DELTA)?"":formatTransCoef(trans);
		return x+y+z+t;

	}

	private static String formatCoef(double c, boolean leading) {
		if (leading) {
			return (deltaComp(Math.abs(c),1,DELTA)?(c>0?"":"-"):String.format("%4.2f",c));
		} else {
			return (deltaComp(Math.abs(c),1,DELTA)?(c>0?"+":"-"):String.format("%+4.2f",c));
		}
	}

	private static String formatTransCoef(double c) {
		if (Math.abs((Math.rint(c)-c))<DELTA) { // this is an integer
			return String.format("%+d",(int)Math.rint(c));
		} else { // it is a fraction
			int num,den;
			int floor = (int)Math.floor(c);
			double decPart = c - floor;
			if (deltaComp(decPart,0.3333333,DELTA)) {
				num=1;den=3;
			} else if (deltaComp(decPart,0.6666667,DELTA)) {
				num=2;den=3;
			} else if (deltaComp(decPart,0.2500000,DELTA)) {
				num=1;den=4;
			} else if (deltaComp(decPart,0.5000000,DELTA)) {
				num=1;den=2;
			} else if (deltaComp(decPart,0.7500000,DELTA)) {
				num=3;den=4;
			} else if (deltaComp(decPart,0.1666667,DELTA)) {
				num=1;den=6;
			} else if (deltaComp(decPart,0.8333333,DELTA)) {
				num=5;den=6;
			} else {
				num=0;den=0; // this in an error
			}
			num = floor*den+num;
			return String.format("%+d/%d", num,den);
			//return String.format("%+4.2f",c);
		}
	}

	protected static boolean deltaComp(double d1, double d2, double delta) {
		return Math.abs(d1-d2)<delta;
	}

	public BravaisLattice getBravLattice() {
		return bravLattice;
	}

	public boolean isEnantiomorphic() {
		Matcher m = nonEnantPat.matcher(shortSymbol);
		if (m.find()) {
			return false;
		}
		return true;
	}

	/**
	 * Given a rotation matrix calculates the rotation axis and angle for it.
	 * The angle is calculated from the trace, the axis from the eigenvalue
	 * decomposition.
	 * If given matrix is improper rotation or identity matrix then
	 * axis (0,0,0) and angle 0 are returned.
	 * @param m
	 * @return
	 * @throws IllegalArgumentException if given matrix is not a rotation matrix (determinant not 1 or -1)
	 */
	public static AxisAngle4d getRotAxisAndAngle(Matrix3d m) {
		double determinant = m.determinant();

		if (!(Math.abs(determinant)-1.0<DELTA)) throw new IllegalArgumentException("Given matrix is not a rotation matrix");

		AxisAngle4d axisAndAngle = new AxisAngle4d(new Vector3d(0,0,0),0);

		double[] d = {m.m00,m.m10,m.m20,
				m.m01,m.m11,m.m21,
				m.m02,m.m12,m.m22};

		Matrix r = new Matrix(d,3);

		if (!deltaComp(r.det(), 1.0, DELTA)) {
			// improper rotation: we return axis 0,0,0 and angle 0
			return axisAndAngle;
		}

		EigenvalueDecomposition evd = new EigenvalueDecomposition(r);

		Matrix eval = evd.getD();
		if (deltaComp(eval.get(0, 0),1.0,DELTA) && deltaComp(eval.get(1, 1),1.0,DELTA) && deltaComp(eval.get(2, 2),1.0,DELTA)) {
			// the rotation is an identity: we return axis 0,0,0 and angle 0
			return axisAndAngle;
		}
		int indexOfEv1;
		for (indexOfEv1=0;indexOfEv1<3;indexOfEv1++) {
			if (deltaComp(eval.get(indexOfEv1, indexOfEv1),1,DELTA)) break;
		}
		Matrix evec = evd.getV();
		axisAndAngle.set(new Vector3d(evec.get(0,indexOfEv1), evec.get(1, indexOfEv1), evec.get(2, indexOfEv1)),
				Math.acos((eval.trace()-1.0)/2.0));

		return axisAndAngle;
	}

	/**
	 * Given a transformation matrix containing a rotation returns the type of rotation:
	 * 1 for identity, 2 for 2-fold rotation, 3 for 3-fold rotation, 4 for 4-fold rotation,
	 * 6 for 6-fold rotation,
	 * -1 for inversions, -2 for mirror planes, -3 for 3-fold improper rotation,
	 * -4 for 4-fold improper rotation and -6 for 6-fold improper rotation
	 * @param m
	 * @return
	 */
	public static int getRotAxisType(Matrix4d m) {
		int axisType = 0;

		Matrix3d rot = new Matrix3d(m.m00,m.m01,m.m02,
				m.m10,m.m11,m.m12,
				m.m20,m.m21,m.m22);

		double determinant = rot.determinant();

		if (!deltaComp(determinant,1.0,DELTA) && !deltaComp(determinant, -1.0, DELTA)) {
			throw new IllegalArgumentException("Given matrix does not seem to be a rotation matrix.");
		}

		int trace = (int)(rot.m00+rot.m11+rot.m22);
		if (determinant>0) {
			switch (trace) {
			case 3:
				axisType=1;
				break;
			case -1:
				axisType=2;
				break;
			case 0:
				axisType=3;
				break;
			case 1:
				axisType=4;
				break;
			case 2:
				axisType=6;
				break;
			default:
				throw new RuntimeException("Trace of transform does not correspond to one of the expected types. This is most likely a bug");
			}
		} else {
			switch (trace) {
			case -3:
				axisType=-1;
				break;
			case 1:
				axisType=-2;
				break;
			case 0:
				axisType=-3;
				break;
			case -1:
				axisType=-4;
				break;
			case -2:
				axisType=-6;
				break;
			default:
				throw new RuntimeException("Trace of transform does not correspond to one of the expected types. This is most likely a bug");
			}
		}
		return axisType;
	}

	@Override
	public String toString() {
		return getShortSymbol();
	}

	public String toXML(){


		JAXBContext jaxbContextStringSortedSet = null;

		try {
			jaxbContextStringSortedSet= JAXBContext.newInstance(SpaceGroup.class);
		} catch (JAXBException e){
			logger.error("Error converting to XML",e);
			return null;
		}

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		PrintStream ps = new PrintStream(baos);


		try {

			Marshaller m = jaxbContextStringSortedSet.createMarshaller();

			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

			m.marshal( this, ps);


		} catch (JAXBException e){
			logger.error("Error converting to XML",e);
		}

		return baos.toString();
	}

	//@XmlElementWrapper(name="transfAlgebraicList", namespace="http://www.biojava.org")
	//@XmlElement
	@XmlJavaTypeAdapter(TransfAlgebraicAdapter.class)
	public List<String> getTransfAlgebraic() {
		return transfAlgebraic;
	}

	public void setTransfAlgebraic(List<String> transfAlgebraic) {
		//System.out.println("setting transfAlgebraic " + transfAlgebraic);
		if ( transformations == null || transformations.size() == 0)
			transformations = new ArrayList<Matrix4d>(transfAlgebraic.size());

		if ( this.transfAlgebraic == null || this.transfAlgebraic.size() == 0)
			this.transfAlgebraic = new ArrayList<String>(transfAlgebraic.size());

		for ( String transf : transfAlgebraic){
			addTransformation(transf);
		}
	}


	public int[] getAxisTypes() {
		return axisTypes;
	}

	public void setAxisTypes(int[] axisTypes) {
		this.axisTypes = axisTypes;
	}

	public static long getSerialversionuid() {
		return serialVersionUID;
	}

	public static Pattern getSplitpat1() {
		return splitPat1;
	}

	public static Pattern getSplitpat2() {
		return splitPat2;
	}

	public static Pattern getCoordpat() {
		return coordPat;
	}

	public static Pattern getTranscoefpat() {
		return transCoefPat;
	}

	public static Pattern getNonenantpat() {
		return nonEnantPat;
	}

	public static double getDelta() {
		return DELTA;
	}

	public void setId(int id) {
		this.id = id;
	}

	public void setMultiplicity(int multiplicity) {

		this.multiplicity = multiplicity;
	}

	public void setPrimitiveMultiplicity(int primitiveMultiplicity) {
		this.primitiveMultiplicity = primitiveMultiplicity;
	}

	public void setShortSymbol(String shortSymbol) {
		this.shortSymbol = shortSymbol;
	}

	public void setAltShortSymbol(String altShortSymbol) {
		this.altShortSymbol = altShortSymbol;
	}


	public void setBravLattice(BravaisLattice bravLattice) {
		this.bravLattice = bravLattice;
	}

}

