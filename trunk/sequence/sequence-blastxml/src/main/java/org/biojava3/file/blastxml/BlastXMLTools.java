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
package org.biojava3.file.blastxml;

import java.io.InputStream;
import java.io.OutputStream;
import java.io.File;
import java.io.Reader;
import java.io.Writer;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

/**
 * Marshals and unmarshals BlastXML output using JAXB's auto-generation 
 * facilities. The DTD it understands can be found in src/main/resources
 * with the name NCBI_BlastOutput.dtd.
 * This module does not comply with the standard BioJava3 way of parsing
 * files. This may change in future, but for now the work of mapping the
 * JAXB auto-generation code into fixed Java interfaces for use within the
 * automated BioJava3 systems would be too much, and would be difficult to
 * maintain to keep up with changes to the DTD over time.
 * @author Mark Schreiber
 * @author Richard Holland
 * @since 3.0
 */
public class BlastXMLTools {

    /**
     * Reads BlastXML output from a source.
     * @param file the source.
     * @return the BlastOutput auto-generated JAXB object representing the
     * BlastXML report. There are no JavaDocs for this object, but it follows
     * the structure of the NCBI_BlastOutput DTD file exactly.
     * @throws JAXBException if parsing failed.
     */
    public BlastOutput unmarshall(File file) throws JAXBException {
        JAXBContext jaxbCtx = JAXBContext.newInstance("org.biojava3.parser.blastxml");
        Unmarshaller unmarshaller = jaxbCtx.createUnmarshaller();
        return (BlastOutput) unmarshaller.unmarshal(file); //NOI18N
    }

    /**
     * Reads BlastXML output from a source.
     * @param is the source.
     * @return the BlastOutput auto-generated JAXB object representing the
     * BlastXML report. There are no JavaDocs for this object, but it follows
     * the structure of the NCBI_BlastOutput DTD file exactly.
     * @throws JAXBException if parsing failed.
     */
    public BlastOutput unmarshall(InputStream is) throws JAXBException {
        JAXBContext jaxbCtx = JAXBContext.newInstance("org.biojava3.parser.blastxml");
        Unmarshaller unmarshaller = jaxbCtx.createUnmarshaller();
        return (BlastOutput) unmarshaller.unmarshal(is); //NOI18N
    }

    /**
     * Reads BlastXML output from a source.
     * @param rdr the source.
     * @return the BlastOutput auto-generated JAXB object representing the
     * BlastXML report. There are no JavaDocs for this object, but it follows
     * the structure of the NCBI_BlastOutput DTD file exactly.
     * @throws JAXBException if parsing failed.
     */
    public BlastOutput unmarshall(Reader rdr) throws JAXBException {
        JAXBContext jaxbCtx = JAXBContext.newInstance("org.biojava3.parser.blastxml");
        Unmarshaller unmarshaller = jaxbCtx.createUnmarshaller();
        return (BlastOutput) unmarshaller.unmarshal(rdr); //NOI18N
    }

    /**
     * Writes BlastXML output to a target.
     * @param blastOutput the BlastOutput auto-generated JAXB object representing the
     * BlastXML report. There are no JavaDocs for this object, but it follows
     * the structure of the NCBI_BlastOutput DTD file exactly.
     * @param file the target.
     * @throws JAXBException if writing failed.
     */
    public void marshall(BlastOutput blastOutput, File file) throws JAXBException {
        JAXBContext jaxbCtx = JAXBContext.newInstance("org.biojava3.parser.blastxml");
        Marshaller marshaller = jaxbCtx.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_ENCODING, "UTF-8"); //NOI18N
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
        marshaller.marshal(blastOutput, file);
    }

    /**
     * Writes BlastXML output to a target.
     * @param blastOutput the BlastOutput auto-generated JAXB object representing the
     * BlastXML report. There are no JavaDocs for this object, but it follows
     * the structure of the NCBI_BlastOutput DTD file exactly.
     * @param rdr the target.
     * @throws JAXBException if writing failed.
     */
    public void marshall(BlastOutput blastOutput, OutputStream os) throws JAXBException {
        JAXBContext jaxbCtx = JAXBContext.newInstance("org.biojava3.parser.blastxml");
        Marshaller marshaller = jaxbCtx.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_ENCODING, "UTF-8"); //NOI18N
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
        marshaller.marshal(blastOutput, os);
    }

    /**
     * Writes BlastXML output to a target.
     * @param blastOutput the BlastOutput auto-generated JAXB object representing the
     * BlastXML report. There are no JavaDocs for this object, but it follows
     * the structure of the NCBI_BlastOutput DTD file exactly.
     * @param wtr the target.
     * @throws JAXBException if writing failed.
     */
    public void marshall(BlastOutput blastOutput, Writer wtr) throws JAXBException {
        JAXBContext jaxbCtx = JAXBContext.newInstance("org.biojava3.parser.blastxml");
        Marshaller marshaller = jaxbCtx.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_ENCODING, "UTF-8"); //NOI18N
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
        marshaller.marshal(blastOutput, wtr);
    }
}
