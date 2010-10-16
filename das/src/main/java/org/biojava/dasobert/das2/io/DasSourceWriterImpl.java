/*
 *                  BioJava development code
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
 * Created on Feb 6, 2006
 *
 */
package org.biojava.dasobert.das2.io;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.dasobert.das.Capabilities;
import org.biojava.dasobert.das.DasTimeFormat;
import org.biojava.dasobert.das2.Das2Capability;
import org.biojava.dasobert.das2.Das2CapabilityImpl;
import org.biojava.dasobert.das2.Das2Source;
import org.biojava.dasobert.dasregistry.Das1Source;
import org.biojava.dasobert.dasregistry.DasCoordinateSystem;
import org.biojava.dasobert.dasregistry.DasSource;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;

public class DasSourceWriterImpl implements DasSourceWriter {

	public static final String COORDSYSURI = "http://www.dasregistry.org/dasregistry/coordsys/";

	public DasSourceWriterImpl() {
		super();

	}

	public void writeCoordinateSystem(XMLWriter xw, DasCoordinateSystem co)
			throws IOException {
		xw.openTag("COORDINATES");
		String uri = co.getUniqueId();
		if (!uri.startsWith(COORDSYSURI))
			uri = COORDSYSURI + uri;

		xw.attribute("uri", uri);

		int taxid = co.getNCBITaxId();
		if (taxid != 0) {
			xw.attribute("taxid", taxid + "");
		}

		xw.attribute("source", co.getCategory());
		xw.attribute("authority", co.getName());
		xw.attribute("test_range", co.getTestCode());
		// TODO: get version from name;
		String version = co.getVersion();
		if ((version != null) && (!version.equals("")))
			xw.attribute("version", version);
		xw.print(co.toString());
		xw.closeTag("COORDINATES");
	}

	public void writeDasSource(XMLWriter xw, DasSource source)
			throws IOException {
		xw.openTag("SOURCE");
		// System.out.println("DasSourceWriterImpl:  writing new source");

		xw.attribute("uri", source.getId());
		xw.attribute("title", source.getNickname());
		String helperurl = source.getHelperurl();
		if (helperurl == null)
			helperurl = "";
		if (!helperurl.equals(""))
			xw.attribute("doc_href", source.getHelperurl());
		String desc = source.getDescription();
		desc = desc.replaceAll("\n", " ");
		desc = desc.replaceAll("\r", " ");
		xw.attribute("description", desc);

		xw.openTag("MAINTAINER");
		xw.attribute("email", source.getAdminemail());
		xw.closeTag("MAINTAINER");
		// System.out.println("before version");
		xw.openTag("VERSION");
		xw.attribute("uri", source.getId());

		Date d = source.getRegisterDate();
		if (d == null)
			d = new Date();
		xw.attribute("created", DasTimeFormat.toDASString(d));

		// System.out.println("before coords");
		List<DasCoordinateSystem> coords = source.getCoordinateSystem();

		for (int i = 0; i < coords.size(); i++) {
			DasCoordinateSystem co = coords.get(i);
			writeCoordinateSystem(xw, co);

		}
		// System.out.println("before das specific part");
		if (source instanceof Das2Source) {
			// System.out.println("das2source");
			Das2Source s = (Das2Source) source;
			Das2Capability[] caps = s.getDas2Capabilities();

			for (int i = 0; i < caps.length; i++) {
				Das2Capability cap = caps[i];
				// System.out.println("DasSourceWriterImpl: capability: " +
				// cap);
				xw.openTag("CAPABILITY");
				xw.attribute("type", cap.getCapability());
				xw.attribute("query_uri", cap.getQueryUri());
				xw.closeTag("CAPABILITY");
			}

		} else if (source instanceof Das1Source) {
			// System.out.println("das1source");
			List<Capabilities> capabilities = source.getCapabilities();
			for (int i = 0; i < capabilities.size(); i++) {
				Capabilities c = capabilities.get(i);
				xw.openTag("CAPABILITY");
				xw.attribute("type", Das2CapabilityImpl.DAS1_CAPABILITY_PREFIX
						+ c);
				xw.attribute("query_uri", source.getUrl() + c);
				xw.closeTag("CAPABILITY");
			}

		}

		List<String> labels = source.getLabels();
		if (labels != null) {

			for (int i = 0; i < labels.size(); i++) {

				xw.openTag("PROP");
				xw.attribute("name", Das2SourceHandler.LABELPROPERTY);
				xw.attribute("value", labels.get(i));
				xw.closeTag("PROP");
			}

		}

		xw.openTag("PROP");
		xw.attribute("name", "leaseTime");

		xw.attribute("value", DasTimeFormat.toDASString(source.getLeaseDate()));
		xw.closeTag("PROP");

		Map<String, String> props = source.getProperties();

		Set<String> keys = props.keySet();
		for (String key : keys) {

			String prop = props.get(key);
			xw.openTag("PROP");
			xw.attribute("name", key);
			xw.attribute("value", prop);
			xw.closeTag("PROP");

		}

		// show validated capabilities
		if (source.getValidCapabilities() != null) {
			List<String> validCapabilities = source.getValidCapabilities();

			// System.out.println("valid caps length="+validCapabilities.length);
			for (int i = 0; i < validCapabilities.size(); i++) {
				xw.openTag("PROP");
				xw.attribute("name", "valid");
				xw.attribute("value", validCapabilities.get(i));

				xw.closeTag("PROP");
			}
		}

		xw.closeTag("VERSION");

		xw.closeTag("SOURCE");
	}

	public void writeDasSource(OutputStream stream, DasSource source)
			throws IOException {
		// System.out.println(source.getNickname());

		PrintWriter pw = new PrintWriter(stream);
		PrettyXMLWriter xw = new PrettyXMLWriter(pw);
		writeDasSource(xw, source);

	}

}
