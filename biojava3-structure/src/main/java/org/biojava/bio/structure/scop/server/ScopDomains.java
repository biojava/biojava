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
package org.biojava.bio.structure.scop.server;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlRootElement;

import org.biojava.bio.structure.scop.ScopDomain;

@XmlRootElement(name = "ScopDomains", namespace ="http://source.rcsb.org")
@XmlAccessorType(XmlAccessType.PUBLIC_MEMBER)
public class ScopDomains implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 7693404355005856746L;

	List<ScopDomain> domains ;
	
	static JAXBContext jaxbContext;
	static {
		try {
			jaxbContext= JAXBContext.newInstance(ScopDomains.class);
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	
	public void setScopDomain(List<ScopDomain> domains) {
		this.domains = domains;
		
	}

	public List<ScopDomain> getScopDomain() {
		return domains;
	}
	
	public  String toXML(){

		ByteArrayOutputStream baos = new ByteArrayOutputStream();

		PrintStream ps = new PrintStream(baos);

		try {

			Marshaller m = jaxbContext.createMarshaller();

			m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

			m.marshal( this, ps);
			

		} catch (Exception e){
			e.printStackTrace();
		}

		return baos.toString();

	}

	public static ScopDomains fromXML(String xml){

		ScopDomains job = null;

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();

			ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

			job = (ScopDomains) un.unmarshal(bais);

		} catch (Exception e){
			e.printStackTrace();
		}

		return job;
	}
	



}
