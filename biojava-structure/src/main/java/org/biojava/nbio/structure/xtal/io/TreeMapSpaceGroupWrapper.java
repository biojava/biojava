/**
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
 * Created on Aug 31, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.xtal.io;

import org.biojava.nbio.structure.xtal.SpaceGroup;

import jakarta.xml.bind.JAXBContext;
import jakarta.xml.bind.JAXBException;
import jakarta.xml.bind.Marshaller;
import jakarta.xml.bind.Unmarshaller;
import jakarta.xml.bind.annotation.XmlAccessType;
import jakarta.xml.bind.annotation.XmlAccessorType;
import jakarta.xml.bind.annotation.XmlRootElement;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.TreeMap;



@XmlRootElement(name = "TreeSetSpaceGroupWrapper", namespace ="http://www.biojava.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class TreeMapSpaceGroupWrapper implements Serializable{


	private static final long serialVersionUID = 4193799052494327416L;

	private TreeMap<Integer,SpaceGroup> data;


	public TreeMapSpaceGroupWrapper(){
		data = new TreeMap<Integer,SpaceGroup>();
	}

	public TreeMap<Integer,SpaceGroup> getData() {
		return data;
	}

	public void setData(TreeMap<Integer,SpaceGroup> data) {
		this.data = data;
	}

	public  String toXML() throws JAXBException{

		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		PrintStream ps = new PrintStream(baos);

		JAXBContext jaxbContext = JAXBContext.newInstance(TreeMapSpaceGroupWrapper.class);

		Marshaller m = jaxbContext.createMarshaller();

		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

		m.marshal( this, ps);

		JAXBContext.newInstance(TreeMapSpaceGroupWrapper.class);

		return baos.toString();

	}

	public static TreeMapSpaceGroupWrapper fromXML(String xml) throws JAXBException{

		TreeMapSpaceGroupWrapper job = null;

		JAXBContext jaxbContext = JAXBContext.newInstance(TreeMapSpaceGroupWrapper.class);

		Unmarshaller un = jaxbContext.createUnmarshaller();

		ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

		job = (TreeMapSpaceGroupWrapper) un.unmarshal(bais);

		return job;
	}


}
