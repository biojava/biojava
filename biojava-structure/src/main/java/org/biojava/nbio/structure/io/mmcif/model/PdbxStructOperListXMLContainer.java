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
package org.biojava.nbio.structure.io.mmcif.model;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.List;

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
