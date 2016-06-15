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

import javax.xml.bind.annotation.adapters.XmlAdapter;
import java.util.Map;
import java.util.TreeMap;

class SpaceGroupMapAdapter extends XmlAdapter<SpaceGroupMapElements[], Map<Integer, SpaceGroup>> {
	@Override
	public SpaceGroupMapElements[] marshal(Map<Integer, SpaceGroup> arg0) throws Exception {
		SpaceGroupMapElements[] mapElements = new SpaceGroupMapElements[arg0.size()];
		int i = 0;
		for (Map.Entry<Integer, SpaceGroup> entry : arg0.entrySet())
			mapElements[i++] = new SpaceGroupMapElements(entry.getKey(), entry.getValue());

		return mapElements;
	}

	@Override
	public Map<Integer, SpaceGroup> unmarshal(SpaceGroupMapElements[] arg0) throws Exception {
		Map<Integer, SpaceGroup> r = new TreeMap<Integer, SpaceGroup>();
		for (SpaceGroupMapElements mapelement : arg0)
			r.put(mapelement.key, mapelement.value);
		return r;
	}
}
