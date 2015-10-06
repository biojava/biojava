package org.biojava.nbio.structure.secstruc;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupIterator;
import org.biojava.nbio.structure.Structure;

/**
 * This class contains methods for obtaining and converting secondary
 * structure information from BioJava Structures.
 * <p>
 * As some examples, a List of all the SSE of a Structure can be obtained, 
 * or the SS of a Structure can be predicted and set to the Groups of the 
 * Structure.
 * 
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class SecStrucTools {

	public static List<SecStrucElement> getSSE(Structure s){
		
		List<SecStrucElement> listSSE = new ArrayList<SecStrucElement>();
		
		GroupIterator iter = new GroupIterator(s);
		while (iter.hasNext()){
			Group g = iter.next();

			if (g instanceof AminoAcid){
				AminoAcid aa = (AminoAcid) g;
				//TODO
				
			}
		}
		
		return listSSE;
	}
}
