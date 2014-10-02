package demo;

import java.io.IOException;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Chain;
//import org.biojava.bio.structure.Group;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.align.util.AtomCache;
import org.biojava.bio.structure.contact.AtomContact;
import org.biojava.bio.structure.contact.AtomContactSet;
//import org.biojava.bio.structure.contact.GroupContact;
import org.biojava.bio.structure.contact.GroupContactSet;
//import org.biojava.bio.structure.io.FileParsingParameters;
import org.biojava3.structure.StructureIO;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DemoContacts {

	private static final Logger logger = LoggerFactory.getLogger(DemoContacts.class);

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

		logger.info("Contacting residues (on CA atoms)");

		for (AtomContact contact:contacts) {
			Atom atom1 = contact.getPair().getFirst();
			Atom atom2 = contact.getPair().getSecond();

			logger.info("{}-{} {}-{} : {}",
					atom1.getGroup().getResidueNumber(),
					atom1.getGroup().getPDBName(),
					atom2.getGroup().getResidueNumber(),
					atom2.getGroup().getPDBName(),
					contact.getDistance());
		}

		logger.info("Total number of atom contacts: {}", contacts.size());

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
		logger.info("Total number of residue contacts: {}", groupContacts.size());

		
		contacts = StructureTools.getAtomsInContact(structure.getChain(0),structure.getChain(1),5.5, false);
		
		logger.info("Contacting residues between 2 first chains (all non-H non-hetatoms)");
		
		for (AtomContact contact:contacts) {
			Atom atom1 = contact.getPair().getFirst();
			Atom atom2 = contact.getPair().getSecond();
			
			logger.info("{}:{}-{}-{} || {}:{}-{}-{} : {}",
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
		
		logger.info("Total number of atom contacts: {}", contacts.size());
		
		groupContacts = new GroupContactSet(contacts);
		logger.info("Total number of residue contacts: {}", groupContacts.size());
	}	
}