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
 * Created on Aug 30, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.nbio.structure.scop.server;

import org.biojava.nbio.structure.scop.ScopDescription;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.List;



@XmlRootElement(name = "ScopDescriptions", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ScopDescriptions implements Serializable{


	private static final long serialVersionUID = 4924350548761431852L;

	static JAXBContext jaxbContext;
	static {
		try {
			jaxbContext= JAXBContext.newInstance(ScopDescriptions.class);
		} catch (Exception e){
			throw new RuntimeException("Could not initialize JAXB context for " + ScopDescriptions.class, e);
		}
	}


	List<ScopDescription> scopDescriptions;

	public List<ScopDescription> getScopDescription() {
		return scopDescriptions;
	}

	public void setScopDescription(List<ScopDescription> descriptions) {
		this.scopDescriptions = descriptions;
	}

	public  String toXML(){

		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		PrintStream ps = new PrintStream(baos);

		try {

			Marshaller m = jaxbContext.createMarshaller();

			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

			m.marshal( this, ps);


		} catch (Exception e){
			throw new RuntimeException("Could not convert " + getClass() + " to XML", e);
		}

		return baos.toString();

	}

	public static ScopDescriptions fromXML(String xml){

		ScopDescriptions job = null;

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();

			ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

			job = (ScopDescriptions) un.unmarshal(bais);

		} catch (Exception e){
			throw new RuntimeException("Could not parse " + ScopDescriptions.class + " from XML", e);
		}

		return job;
	}


}
