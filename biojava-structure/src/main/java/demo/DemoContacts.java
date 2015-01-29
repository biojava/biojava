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
package demo;

import org.biojava.bio.structure.*;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.contact.AtomContact;
import org.biojava.bio.structure.contact.AtomContactSet;
import org.biojava.bio.structure.contact.GroupContactSet;
import org.biojava3.structure.StructureIO;

import java.io.IOException;

//import org.biojava.bio.structure.Group;
//import org.biojava.bio.structure.contact.GroupContact;
//import org.biojava.bio.structure.io.FileParsingParameters;

public class DemoContacts {

	
	public static void main(String[] args) throws IOException, StructureException {

		String pdbCode = "1smt";
		
		demoContacts(pdbCode);
	}
	
	private static void demoContacts(String pdbCode) throws IOException, StructureException {
		
		AtomCache cache = new AtomCache();
		cache.setUseMmCif(true);
		
		StructureIO.setAtomCache(cache); 
		
		Structure structure = StructureIO.getStructure(pdbCode);
			
		Chain chain = structure.getChainByPDB("A");
		
			
		String[] atoms = {"CA"};
		AtomContactSet contacts = StructureTools.getAtomsInContact(chain, atoms, 8.0);

		System.out.println("Contacting residues (on CA atoms)");

		for (AtomContact contact:contacts) {
			Atom atom1 = contact.getPair().getFirst();
			Atom atom2 = contact.getPair().getSecond();

			System.out.printf(" %3s-%3s %3s-%3s : %5.2f\n",
					atom1.getGroup().getResidueNumber(),
					atom1.getGroup().getPDBName(),
					atom2.getGroup().getResidueNumber(),
					atom2.getGroup().getPDBName(),
					contact.getDistance());
		}

		System.out.println("Total number of atom contacts: "+contacts.size());

		GroupContactSet groupContacts = new GroupContactSet(contacts);
//		for (GroupContact groupContact:groupContacts) {
//			Group g1 = groupContact.getPair().getFirst();
//			Group g2 = groupContact.getPair().getSecond();
//
//			System.out.printf(" %3s-%3s %3s-%3s : %5.2f\n",
//					g1.getResidueNumber(),
//					g1.getPDBName(),
//					g2.getResidueNumber(),
//					g2.getPDBName(),
//					groupContact.getMinDistance());
//		}
		System.out.println("Total number of residue contacts: "+groupContacts.size());

		
		contacts = StructureTools.getAtomsInContact(structure.getChain(0),structure.getChain(1),5.5, false);
		
		System.out.println("Contacting residues between 2 first chains (all non-H non-hetatoms)");
		
		for (AtomContact contact:contacts) {
			Atom atom1 = contact.getPair().getFirst();
			Atom atom2 = contact.getPair().getSecond();
			
			System.out.printf(" %3s:%1s-%3s-%3s || %3s:%1s-%3s-%3s : %5.2f\n",
					atom1.getGroup().getResidueNumber(),
					atom1.getGroup().getChainId(),
					atom1.getGroup().getPDBName(),
					atom1.getName(),
					atom2.getGroup().getResidueNumber(),
					atom2.getGroup().getChainId(),
					atom2.getGroup().getPDBName(),
					atom2.getName(),
					contact.getDistance());
		}
		
		System.out.println("Total number of atom contacts: "+contacts.size());
		
		groupContacts = new GroupContactSet(contacts);
		System.out.println("Total number of residue contacts: "+groupContacts.size());

	}

	
}
