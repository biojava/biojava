package org.biojava.bio.structure.io.mmcif.model;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name="PdbxStructAssemblyXMLContainer")
public class PdbxStructAssemblyXMLContainer {

	private List<PdbxStructAssembly> data ;
	
	static JAXBContext jaxbContext;
	static {
		try {
			jaxbContext= JAXBContext.newInstance(PdbxStructAssemblyXMLContainer.class);
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	@XmlElementWrapper
	public List<PdbxStructAssembly> getPdbxStructAssemblies(){
		return data;
		
	}
	
	public void setPdbxStructAssemblies(List<PdbxStructAssembly> d){
		data = d;
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

	public static PdbxStructAssemblyXMLContainer fromXML(String xml){

		PdbxStructAssemblyXMLContainer job = null;

		try {

			Unmarshaller un = jaxbContext.createUnmarshaller();

			ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

			job = (PdbxStructAssemblyXMLContainer) un.unmarshal(bais);

		} catch (Exception e){
			e.printStackTrace();
		}

		return job;
	}
	
}
