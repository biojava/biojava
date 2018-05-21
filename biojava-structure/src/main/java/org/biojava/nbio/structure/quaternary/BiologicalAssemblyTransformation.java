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
package org.biojava.nbio.structure.quaternary;

import org.biojava.nbio.structure.xtal.CrystalCell;
import org.biojava.nbio.structure.xtal.CrystalTransform;
import org.biojava.nbio.core.util.PrettyXMLWriter;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * The transformation needed for generation of biological assemblies
 * from the contents of a PDB/mmCIF file. It contains both the actual
 * transformation (rotation+translation) and the chain identifier to
 * which it should be applied.
 *
 * @author Peter Rose
 * @author Andreas Prlic
 * @author rickb
 * @author Jose Duarte
 * @see CrystalTransform
 */
public class BiologicalAssemblyTransformation implements Cloneable, Comparable<BiologicalAssemblyTransformation>, Serializable {

	private static final long serialVersionUID = -6388503076022480391L;

	private String id;
	private String chainId;
	private Matrix4d transformation;

	/**
	 * Default Constructor
	 */
	public BiologicalAssemblyTransformation() {
		// we initialize to identity so that setting rotation and translation work properly
		transformation = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);
	}

	/**
	 * Copy Constructor
	 *
	 * @param src
	 */
	public BiologicalAssemblyTransformation(final BiologicalAssemblyTransformation src)
	{
		this.transformation = new Matrix4d(src.transformation);
		this.id = src.getId();
		this.chainId =  src.getChainId();
	}


	/**
	 * Sets the identifier for this biological assembly transformation. This is usually
	 * the model number used in the biological assembly files.
	 * @param id
	 */
	public void setId(String id) {
		this.id = id;
	}

	/**
	 * Returns the identifier for this biological assembly transformation.
	 * @return biological assembly transformation identifier
	 */
	public String getId() {
		return id;
	}

	/**
	 * Sets the chain identifier (asym id) that this transformation should be applied to.
	 * @param chainId
	 */
	public void setChainId(String chainId) {
		this.chainId = chainId;
	}

	/**
	 * Returns the chain identifier (asym id) that this transformation should be applied to.
	 * @return chain identifier
	 */
	public String getChainId() {
		return this.chainId;
	}

	/**
	 * Sets the transformation using a 4x4 transformation matrix
	 * @param transformation
	 */
	public void setTransformationMatrix(Matrix4d transformation) {
		this.transformation = transformation;
	}

	/**
	 * Return the transformation (both rotational and translational component) as a 4x4 transformation matrix.
	 * The transformation is in orthonormal (cartesian coordinates). If required to be converted to
	 * crystal coordinates then use {@link CrystalCell#transfToCrystal(Matrix4d)}
	 * Note that this is a reference to the variable, thus it remains linked to this object's transformation field.
	 * The user must deep copy it if need changing it.
	 * @return 4x4 transformation matrix
	 */
	public Matrix4d getTransformationMatrix() {
		return transformation;
	}

	public void setRotationMatrix(double[][] m) {
		for (int i=0;i<3;i++) {
			for (int j=0;j<3;j++) {
				this.transformation.setElement(i, j, m[i][j]);
			}
		}
	}

	public void setTranslation(double[] t) {
		for (int i=0;i<3;i++) {
			this.transformation.setElement(i, 3, t[i]);
		}
	}

	/**
	 * Applies the transformation to given point.
	 */
	public void transformPoint(final double[] point) {
		Point3d p = new Point3d(point[0],point[1],point[2]);
		transformation.transform(p);
		point[0] = p.x;
		point[1] = p.y;
		point[2] = p.z;
	}

	/**
	 * Returns the combination (product) of two biological assembly transformations.
	 * @param matrix1
	 * @param matrix2
	 * @return combined transformation
	 */
	public static BiologicalAssemblyTransformation combine(BiologicalAssemblyTransformation matrix1, BiologicalAssemblyTransformation matrix2) {
		Matrix4d transformation = new Matrix4d(matrix1.transformation);
		transformation.mul(matrix2.transformation);
		BiologicalAssemblyTransformation combined = new BiologicalAssemblyTransformation();
		combined.setTransformationMatrix(transformation);
		return combined;
	}
	
	/**
	 * Tells whether this transformation is in identity.
	 * @return
	 */
	public boolean isIdentity() {
		return transformation.epsilonEquals(new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1), 0.00000000001);
	}

	public String toXML() throws IOException{

		StringWriter sw = new StringWriter();
		PrintWriter writer = new PrintWriter(sw);

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(writer));

		toXML(xml);

		xml.close();
		writer.close();
		sw.close();
		return sw.toString();
	}

	public void toXML(PrettyXMLWriter xml) throws IOException{
		xml.openTag("transformation");
		xml.attribute("index",id);

		xml.openTag("matrix");

			for ( int i = 0 ; i<3 ; i++){
				for ( int j = 0 ; j<3 ;j++){
				xml.attribute("m" +  (i+1) + (j+1), String.format("%.8f",transformation.getElement(i,j)));
			}
		}
		xml.closeTag("matrix");

		xml.openTag("shift");
		for ( int i = 0 ; i<3 ; i++) {
			xml.attribute("v"+(i+1),String.format("%.8f", transformation.getElement(i,3)));
		}
		xml.closeTag("shift");

		xml.closeTag("transformation");

	}

	public static BiologicalAssemblyTransformation fromXML(String xml)
			throws SAXException,
			IOException,
			ParserConfigurationException{


		List<BiologicalAssemblyTransformation> transformations = fromMultiXML(xml);

		if ( transformations.size() > 0)
			return transformations.get(0);

		else
			return null;
	}

	public static List<BiologicalAssemblyTransformation> fromMultiXML(String xml) throws ParserConfigurationException, SAXException, IOException{


		List<BiologicalAssemblyTransformation> transformations = new ArrayList<BiologicalAssemblyTransformation>();

		// read the XML of a string and returns a ModelTransformationmatrix
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = factory.newDocumentBuilder();

		InputSource inStream = new InputSource();
		inStream.setCharacterStream(new StringReader(xml));
		Document doc = db.parse(inStream);

		// normalize text representation
		doc.getDocumentElement().normalize();;

		NodeList listOfTransforms = doc.getElementsByTagName("transformation");

		for(int pos=0; pos<listOfTransforms.getLength() ; pos++) {
			Node rootElement       = listOfTransforms.item(pos);

			BiologicalAssemblyTransformation max = new BiologicalAssemblyTransformation();

			max.id = getAttribute(rootElement,"index");
			max.chainId = getAttribute(rootElement,"chainName");

			NodeList listOfChildren = rootElement.getChildNodes();


			for(int i=0; i<listOfChildren.getLength() ; i++)
			{
				// and now the matrix ...
				Node block       = listOfChildren.item(i);

				// we only look at blocks.
				if ( block.getNodeName().equals("matrix"))
					max.setRotationMatrix(getMatrixFromXML(block));

				if ( block.getNodeName().equals("shift"))
					max.setTranslation(getVectorFromXML(block));

			}

			transformations.add(max);
		}

		return transformations;
	}

	private static double[] getVectorFromXML(Node block) {
		double[] d = new double[3];
		for ( int i = 0 ; i<3 ; i++){
			d[i] = Float.parseFloat(getAttribute(block, "v" + (i+1) ));
		}
		return d;
	}

	private static double[][] getMatrixFromXML(Node block) {
		double[][] m  = new double[3][3];
		for ( int i = 0 ; i<3 ; i++){
			for ( int j = 0 ; j<3 ;j++){
				// TODO check is this matrix is populated correctly
				String val = getAttribute(block, "m" + (j+1)+(i+1));
				m[j][i] = Double.parseDouble(val);
			}
		}
		return m;
	}

	private static String getAttribute(Node node, String attr){
		if( ! node.hasAttributes())
			return null;

		NamedNodeMap atts = node.getAttributes();

		if ( atts == null)
			return null;

		Node att = atts.getNamedItem(attr);
		if ( att == null)
			return null;

		String value = att.getTextContent();

		return value;

	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "BiologicalAssemblyTransformation [id=" + id + ", chainId="
				+ chainId + ", rotation=" + rotMatrixToString(transformation) + ", translation="
				+ translVecToString(transformation) + "]";
	}

	public static String rotMatrixToString(Matrix4d m) {
		return String.format("(%5.2f %5.2f %5.2f, %5.2f %5.2f %5.2f, %5.2f %5.2f %5.2f)", m.m00, m.m01, m.m02,  m.m10, m.m11, m.m12,  m.m20, m.m21, m.m22);
	}

	public static String translVecToString(Matrix4d m) {
		return String.format("(%5.2f %5.2f %5.2f)", m.m03, m.m13, m.m23);
	}

    @Override
    public int compareTo(BiologicalAssemblyTransformation other) {
        int comp = this.chainId.compareTo(other.chainId);
        return comp == 0 ? this.id.compareTo(other.id) : comp;
    }
}
