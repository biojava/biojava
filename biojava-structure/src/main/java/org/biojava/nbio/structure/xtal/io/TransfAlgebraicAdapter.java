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
 * Created on Mar 29, 2014
 * Author: andreas
 *
 */

package org.biojava.nbio.structure.xtal.io;

import javax.xml.bind.annotation.adapters.XmlAdapter;
import java.util.ArrayList;
import java.util.List;

public class TransfAlgebraicAdapter extends XmlAdapter<String[], List<String>>{

	@Override
	public String[] marshal(List<String> arg0) throws Exception {
		String[] elements = new String[arg0.size()];
		int i = 0;
		for (String s : arg0)
			elements[i++] = s;
		return elements;
	}

	@Override
	public List<String> unmarshal(String[] arg0) throws Exception {
		List<String> l = new ArrayList<String>();
		for (String s : arg0)
			l.add(s);
		return l;
	}

}
