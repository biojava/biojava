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
package org.biojava.nbio.structure.io;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Bond;
import org.biojava.nbio.structure.Group;


public class GroupToSDF {


	public String getText(Group thisGroup) {
		// Fuction to convert a Group to a strign of  the MDL molnlock
		StringBuilder sb = new StringBuilder();
		sb.append(getHeader(thisGroup));
		sb.append(getCtab(thisGroup));
		return sb.toString();
	}

	private String getCtab(Group thisGroup){
		DecimalFormat df = new DecimalFormat("0.0000");
		// The thre string builders for the three parts of the MDL block
		StringBuilder header = new StringBuilder();
		StringBuilder atomList = new StringBuilder();
		StringBuilder bondOrders = new StringBuilder();
		int numBonds = 0;
		List<Atom> atoms = thisGroup.getAtoms();
		for(Atom a: thisGroup.getAtoms()){
			/// ALL SHOULD BE TO FOUR DECIMAL PLACES WITH A CERTAIN NUMBER OF SPACES BEFOREa
			String spaceX = getSpace(10, df.format(a.getX()));
			String spaceY = getSpace(10, df.format(a.getY()));
			String spaceZ = getSpace(10, df.format(a.getZ()));
			String spaceEle = getSpace(4, a.getElement().toString());
			atomList.append(spaceX+df.format(a.getX())+spaceY+df.format(a.getY())+spaceZ+df.format(a.getZ())+" "+a.getElement().toString()+spaceEle+"0  0  0  0  0  0  0  0  0  0  0  0\n");
			for(Bond b: a.getBonds()){
				Atom otherAtom = b.getOther(a);
				if(atoms.indexOf(otherAtom)>=atoms.indexOf(a)){
					continue;
				}
				if(atoms.indexOf(otherAtom)<0){
					continue;
				}
				if(atoms.indexOf(a)<0){
					continue;
				}
				numBonds++;
				// Deal with the index infromation
				String indexOther = Integer.toString(atoms.indexOf(otherAtom)+1);
				String index = Integer.toString(atoms.indexOf(a)+1);
				String order = Integer.toString(b.getBondOrder());
				String spaceIndOne = getSpace(3, index);
				String spaceIndTwo = getSpace(3, indexOther);
				String spaceOrder = getSpace(3, order);
				bondOrders.append(spaceIndOne).append(index).append(spaceIndTwo).append(indexOther).append(spaceOrder).append(order).append("  0\n");
			}
		}
		StringBuilder outString = new  StringBuilder();
		// Add the first line now
		String spaceNumAtoms = getSpace(3, Integer.toString(thisGroup.getAtoms().size()));
		String spaceNumBonds = getSpace(3, Integer.toString(numBonds));
		header.append("\n");
		header.append(spaceNumAtoms).append(thisGroup.getAtoms().size()).append(spaceNumBonds).append(numBonds).append("  0  0  0  0  0  0  0  0999 V2000\n");
		// Now add the header, atom, bond and charge information togeyher
		outString.append(header.toString());
		outString.append(atomList.toString());
		outString.append(bondOrders.toString());
		outString.append(getCharges(thisGroup));
		// Add the final line and the $$$$ delimiter
		outString.append("M  END\n$$$$");
		// Get the string and return it
		return outString.toString();
	}

	private Object getCharges(Group thisGroup) {
		// Get the array of charges
		List<Number> chargeList = getAtomCharges(thisGroup);
		// Build the string
		StringBuilder outS = new StringBuilder();
		int chargeCount=0;
		// Loop through the charges -> maximum 8 charges per line
		for(int i =0; i<chargeList.size();i++){
			short charge = (Short) chargeList.get(i);
			if(charge!=0){
				if(chargeCount==0){
					outS.append("M  CHG   N");
				}
				outS.append(getSpace(4, Integer.toString(i))+(i+1));
				outS.append(getSpace(4, Short.toString(charge))+charge);
				chargeCount++;
			}
			//
			if(chargeCount==8){
				outS.append("\n");
				outS.replace(0, 10, "M  CHG   "+chargeCount);
				chargeCount=0;
			}

		}
		if(chargeCount==0){
			return "";
		}
		// Now add the charge count
		outS.replace(0, 10, "M  CHG  "+chargeCount);
		outS.append("\n");
		// Now return the string
		return outS.toString();
	}

	private static List<Number> getAtomCharges(Group group) {
		// The list to store the answer
		List<Number> outArr = new ArrayList<Number>();
		// Get the atom charge Information
		for(Atom a: group.getAtoms()){
			outArr.add(a.getCharge());
		}
		return outArr;
	}

	private String getSpace(int inputNum, String format) {
		// Function to get the spaces between numbers
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<(inputNum-format.length());i++){
			sb.append(" ");
		}
		return sb.toString();
	}

	private String getHeader(Group thisGroup) {
		// Make the header info for the start of the block
		StringBuilder sb = new StringBuilder();
		sb.append(thisGroup.getPDBName()).append("\n");
		sb.append("Made by BioJava");
		sb.append("\n");
		return sb.toString();
	}
}
