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
 * Created on Oct 10, 2011
 * Created by Andreas Prlic
 *
 * @since 3.0.2
 */
package org.biojava.bio.structure.align.xml;

import java.io.StringWriter;
import java.util.SortedSet;
import java.util.TreeSet;

import org.biojava.bio.structure.align.client.PdbPair;

public class PdbPairsMessage {
	
	String method ;
	
	SortedSet<PdbPair> pairs;
	
	public PdbPairsMessage(){
	
		method = PdbPairXMLConverter.DEFAULT_METHOD_NAME;
		
		pairs = new TreeSet<PdbPair>();
	
	}
	
	public String getMethod() {
		return method;
	}
	public void setMethod(String method) {
		this.method = method;
	}
	public SortedSet<PdbPair> getPairs() {
		return pairs;
	}
	public void setPairs(SortedSet<PdbPair> pairs) {
		this.pairs = pairs;
	}
	
	
	public String toString(){
		
		StringWriter w = new StringWriter();
		
		w.append("PdbPairsMessage: ");
		w.append("algorithm: ");	
		w.append(method);
		w.append(" pairs: ");
		w.append(pairs.toString());
			
		
		return w.toString();
		
	}
	
	
	
}
