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

@XmlRootElement(name="PdbxStructOperListXMLContainer")
public class PdbxStructOperListXMLContainer {
	
	

		private List<PdbxStructOperList> data ;
		
		static JAXBContext jaxbContext;
		static {
			try {
				jaxbContext= JAXBContext.newInstance(PdbxStructOperList.class);
			} catch (Exception e){
				e.printStackTrace();
			}
		}
		
		@XmlElementWrapper
		public List<PdbxStructOperList> getPdbxStructOperLists(){
			return data;
			
		}
		
		public void setPdbxStructOperLists(List<PdbxStructOperList> d){
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

		public static PdbxStructOperListXMLContainer fromXML(String xml){

			PdbxStructOperListXMLContainer job = null;

			try {

				Unmarshaller un = jaxbContext.createUnmarshaller();

				ByteArrayInputStream bais = new ByteArrayInputStream(xml.getBytes());

				job = (PdbxStructOperListXMLContainer) un.unmarshal(bais);

			} catch (Exception e){
				e.printStackTrace();
			}

			return job;
		}
		
}
