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
package org.biojava.nbio.structure.io.mmcif;

import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.mmcif.chem.ResidueType;
import org.biojava.nbio.structure.io.mmcif.model.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class ChemCompConsumer implements MMcifConsumer {

	private static final Logger logger = LoggerFactory.getLogger(ChemCompConsumer.class);
	
	ChemicalComponentDictionary dictionary;

	String latestChemCompId;
	public ChemCompConsumer(){
		dictionary = new ChemicalComponentDictionary();
	}

	@Override
	public void documentStart() {


	}

	public ChemicalComponentDictionary getDictionary(){
		return dictionary;
	}

	@Override
	public void newChemComp(ChemComp c) {

		if ( c.getId() == null)
			logger.warn("chem comp ID == null " + c);

		latestChemCompId = c.getId();
		dictionary.addChemComp(c);		
		if ( c.getResidueType() == ResidueType.nonPolymer)
			return;

		if ( c.getResidueType() == ResidueType.saccharide)
			return;

		if ( c.getResidueType() == ResidueType.dSaccharide)
			return;

		//if ( c.isStandard())
		//	System.out.println(c);
	}

	@Override
	public void documentEnd() {


	}

	@Override
	public void newAtomSite(AtomSite atom) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newDatabasePDBremark(DatabasePDBremark remark) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newDatabasePDBrev(DatabasePDBrev dbrev) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newDatabasePDBrevRecord(DatabasePdbrevRecord dbrev) {

	}

	@Override
	public void newEntity(Entity entity) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newEntityPolySeq(EntityPolySeq epolseq) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newExptl(Exptl exptl) {
		// TODO Auto-generated method stub

	}
	
	@Override
	public void newCell(Cell cell) {
		// TODO Auto-generated method stub
	}
	
	@Override
	public void newSymmetry(Symmetry symmetry) {
		// TODO Auto-generated method stub 
	}
	
	@Override
	public void newStructNcsOper(StructNcsOper sNcsOper) {
		// TODO Auto-generated method stub
	}

	@Override
	public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newRefine(Refine r) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newStructAsym(StructAsym sasym) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newStructKeywords(StructKeywords kw) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newStructRef(StructRef sref) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newStructRefSeq(StructRefSeq sref) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newStructRefSeqDif(StructRefSeqDif sref) {
		// TODO Auto-generated method stub

	}

	@Override
	public void setStruct(Struct struct) {
		// TODO Auto-generated method stub

	}

	@Override
	public void newGenericData(String category, List<String> loopFields,
			List<String> lineData) {
		//System.out.println("unhandled category: " + category);

	}


	@Override
	public void newAuditAuthor(AuditAuthor aa)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public FileParsingParameters getFileParsingParameters()
	{
		// can be ingored in this case...
		return null;
	}

	@Override
	public void setFileParsingParameters(FileParsingParameters params)
	{
		// TODO Auto-generated method stub

	}

	@Override
	public void newChemCompDescriptor(ChemCompDescriptor ccd) {
		ChemComp cc = dictionary.getChemComp(latestChemCompId);
		cc.getDescriptors().add(ccd);

	}

	@Override
	public void newPdbxStructOperList(PdbxStructOperList structOper) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxStrucAssembly(PdbxStructAssembly strucAssembly) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxStrucAssemblyGen(PdbxStructAssemblyGen strucAssembly) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newChemCompAtom(ChemCompAtom atom) {
		dictionary.getChemComp(latestChemCompId).getAtoms().add(atom);
	}

	@Override
	public void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newChemCompBond(ChemCompBond bond) {
		dictionary.getChemComp(latestChemCompId).getBonds().add(bond);
	}

	@Override
	public void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newEntitySrcGen(EntitySrcGen entitySrcGen) {
		// TODO Auto-generated method stub
		
	}
	@Override
	public void newEntitySrcNat(EntitySrcNat entitySrcNat) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newEntitySrcSyn(EntitySrcSyn entitySrcSyn) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newStructConn(StructConn structConn) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newStructSiteGen(StructSiteGen gen) {
		// TODO
	}

	@Override
	public void newStructSite(StructSite site) {
		// TODO
	}
}

