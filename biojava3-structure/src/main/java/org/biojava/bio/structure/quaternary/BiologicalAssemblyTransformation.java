package org.biojava.bio.structure.quaternary;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.core.util.PrettyXMLWriter;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * @author Peter Rose
 * @author Andreas Prlic
 * @author rickb
 *
 */
public class BiologicalAssemblyTransformation implements Cloneable, Serializable {

	private static final long serialVersionUID = -6388503076022480391L;
	private String id = null;
	private String chainId = null;
	private Matrix rotation = null;
	private double[] translation = new double[3];
	
	/**
	 * Default Constructor
	 */
	public BiologicalAssemblyTransformation() {}

	/**
	 * Copy Constructor
	 * 
	 * @param src
	 */
	public BiologicalAssemblyTransformation(final BiologicalAssemblyTransformation src)
	{
		this.rotation = new Matrix(src.getRotationMatrix().getArrayCopy());
		this.translation [0] = src.translation[0];
		this.translation [1] = src.translation[1];
		this.translation [2] = src.translation[2];
		this.id = src.getId();
		this.chainId =  src.getChainId();
	}
	
	
	/**
	 * Sets the identifier for this biological assembly transformation. This is usually the model number used in the biological assembly files.
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
	 * Sets the chain identified this transformation should be applied to.
	 * @param chainId
	 */
	public void setChainId(String chainId) {
		this.chainId = chainId;
	}
	
	/**
	 * Returns the chain identifier this transformation should be applied to.
	 * @return chain identifier
	 */
	public String getChainId() {
		return this.chainId;
	}

	/**
	 * Set the rotation matrix for this transformation
	 * @param rotation
	 */
	public void setRotationMatrix(Matrix rotation) {
		this.rotation = new Matrix(rotation.getArray());
	}
	
	/**
	 * Return the rotation matrix for this transformation
	 * @return 3x3 rotation matrix
	 */
	public Matrix getRotationMatrix() {
		return this.rotation;
	}
	
	/**
	 * Sets the translational component of this transformation
	 * @param translation
	 */
	public void setTranslation(double[] translation) {
		this.translation [0] = translation[0];
		this.translation [1] = translation[1];
		this.translation [2] = translation[2];
	}
	
	/**
	 * Returns the translational component of this transformation
	 * @return translational component
	 */
	public double[] getTranslation() {
		return this.translation;
	}
	
	/**
	 * Sets the transformation using a 4x4 transformation matrix
	 * @param transformation
	 */
	public void setTransformationMatrix(Matrix transformation) {
		this.rotation = transformation.getMatrix(0, 2, 0, 2);
		this.translation = new double[3];
		this.translation[0] = transformation.get(0, 3);
		this.translation[1] = transformation.get(1, 3);
		this.translation[2] = transformation.get(2, 3);
	}
	
	/**
	 * Returns the rotational and translational component of this transformation as 4x4 transformation matrix.
	 * @return 4X4 rotation matrix
	 */
	public Matrix getTransformationMatrix() {
		Matrix transformation = Matrix.identity(4, 4);
		transformation.setMatrix(0, 2, 0, 2, rotation);
		transformation.set(0, 3, translation[0]);
		transformation.set(1, 3, translation[1]);
		transformation.set(2, 3, translation[2]);
		return transformation;
	}
	
	/**
	 * Applies the transformation to this point.
	 */
	public void transformPoint(final double[] point) {
		double x = point[0] * rotation.get(0, 0) + point[1] * rotation.get(0, 1) + point[2] * rotation.get(0, 2) + translation[0];
		double y = point[0] * rotation.get(1, 0) + point[1] * rotation.get(1, 1) + point[2] * rotation.get(1, 2) + translation[1];
		double z = point[0] * rotation.get(2, 0) + point[1] * rotation.get(2, 1) + point[2] * rotation.get(2, 2) + translation[2];
		point[0] = x;
		point[1] = y;
		point[2] = z;
	}
	
	/**
	 * Returns the combination (product) of two biological assembly transformations.
	 * @param matrix1 4x4 transformation matrix
	 * @param matrix2 4x4 transformation matrix
	 * @return combined transformation
	 */
	public static BiologicalAssemblyTransformation combine(BiologicalAssemblyTransformation matrix1, BiologicalAssemblyTransformation matrix2) {	
		Matrix transformation = matrix1.getTransformationMatrix().times(matrix2.getTransformationMatrix());
		BiologicalAssemblyTransformation combined = new BiologicalAssemblyTransformation();
		combined.setTransformationMatrix(transformation);
		return combined;
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
		Matrix m = getRotationMatrix();
		
			for ( int i = 0 ; i<3 ; i++){
				for ( int j = 0 ; j<3 ;j++){	
				xml.attribute("m" +  (i+1) + (j+1), String.format("%.8f",m.get(i,j)));
			}
		}
		xml.closeTag("matrix");

		double[] shift = getTranslation();
		xml.openTag("shift");
		for ( int i = 0 ; i<3 ; i++) {
			xml.attribute("v"+(i+1),String.format("%.8f", shift[i]));
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
			max.chainId = getAttribute(rootElement,"chainId");

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

	private static Matrix getMatrixFromXML(Node block) {
		Matrix m  = new Matrix(3,3);
		for ( int i = 0 ; i<3 ; i++){
			for ( int j = 0 ; j<3 ;j++){
				// TODO check is this matrix is populated correctly
				String val = getAttribute(block, "m" + (j+1)+(i+1));			
				m.set(j,i, Float.parseFloat(val));
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
				+ chainId + ", rotation=" + rotation + ", translation="
				+ Arrays.toString(translation) + "]";
	}
	
	

}