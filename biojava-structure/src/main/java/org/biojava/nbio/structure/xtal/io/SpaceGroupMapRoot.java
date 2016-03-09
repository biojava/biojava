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
 * created at 27 Mar 2014
 * Author: ap3
 */

package org.biojava.nbio.structure.xtal.io;

import org.biojava.nbio.structure.xtal.SpaceGroup;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.TreeMap;

@XmlRootElement(name="SpaceGroupMapRoot", namespace ="http://www.biojava.org")

public class SpaceGroupMapRoot {

	private TreeMap<Integer, SpaceGroup> mapProperty;

	public SpaceGroupMapRoot() {
		mapProperty = new TreeMap<Integer, SpaceGroup>();
	}

	@XmlJavaTypeAdapter(SpaceGroupMapAdapter.class)
	public TreeMap<Integer, SpaceGroup> getMapProperty() {
		return mapProperty;
	}

	public void setMapProperty(TreeMap<Integer, SpaceGroup> map) {
		this.mapProperty = map;
	}


	public  String toXML() throws JAXBException {

		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		PrintStream ps = new PrintStream(baos);

		JAXBContext jaxbContext = JAXBContext.newInstance(SpaceGroupMapRoot.class);

		Marshaller xmlConverter = jaxbContext.createMarshaller();
		xmlConverter.setProperty("jaxb.formatted.output",true);
		xmlConverter.marshal(this,ps);

		return baos.toString();
	}

	public static SpaceGroupMapRoot fromXML(String xml) throws JAXBException{
		SpaceGroupMapRoot job = null;

		JAXBContext jaxbContext = JAXBContext.newInstance(SpaceGroupMapRoot.class);

		Unmarshaller un = jaxbContext.createUnmarshaller();

		ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

		job = (SpaceGroupMapRoot) un.unmarshal(bais);


		return job;
	}
}
