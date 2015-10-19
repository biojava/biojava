package org.biojava.nbio.structure.secstruc;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupIterator;
import org.biojava.nbio.structure.ResidueNumber;
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

	/**
	 * Obtain a List of SSE of the Structure, with the proper IDs
	 * and residue ranges.
	 * 
	 * @param s Structure with SS assignments
	 * @return List of SecStrucElement objects
	 */
	public static List<SecStrucElement> getSSE(Structure s){
		
		List<SecStrucElement> listSSE = new ArrayList<SecStrucElement>();
		GroupIterator iter = new GroupIterator(s);
		
		//SecStruc information - initialize
		SecStrucType type = SecStrucType.coil;
		ResidueNumber previous = new ResidueNumber();
		ResidueNumber start = new ResidueNumber();
		String chainId = "";
		int count = 0; //counts the number of residues in SSE
		
		//Create a map for the IDs of the SSE in the structure
		Map<SecStrucType,Integer> ids = new TreeMap<SecStrucType,Integer>();
		for (SecStrucType t : SecStrucType.values()) ids.put(t, 1);
		
		while (iter.hasNext()){
			Group g = iter.next();
			if (g.hasAminoAtoms()){
				SecStrucInfo ss = (SecStrucInfo) g.getProperty(Group.SEC_STRUC);
				if (ss == null) continue;
				if (count > 0){
					if (ss.type == type && chainId == g.getChainId()) {
						previous = g.getResidueNumber();
						count++;
						continue;
					}
					else {
						//Save the current
						SecStrucElement sse = new SecStrucElement(
								type, start, previous, 
								count, ids.get(type), chainId);
						listSSE.add(sse);
						ids.put(type, ids.get(type)+1);
						count = 0;
						
						//Initialize a new SSE one
						if (ss.type != SecStrucType.coil){
							type = ss.type;
							start = g.getResidueNumber();
							previous = start;
							chainId = g.getChainId();
							count = 1;
						}
					}
				} else {
					if (ss.type != SecStrucType.coil){
						type = ss.type;
						start = g.getResidueNumber();
						previous = start;
						chainId = g.getChainId();
						count = 1;
					}
				}
			}
		}
		
		return listSSE;
	}
}
