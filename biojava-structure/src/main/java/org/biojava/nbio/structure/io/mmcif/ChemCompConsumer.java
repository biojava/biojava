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


	}

	@Override
	public void newDatabasePDBremark(DatabasePDBremark remark) {


	}

	@Override
	public void newDatabasePDBrev(DatabasePDBrev dbrev) {


	}

	@Override
	public void newDatabasePDBrevRecord(DatabasePdbrevRecord dbrev) {

	}

	@Override
	public void newEntity(Entity entity) {


	}

	@Override
	public void newEntityPolySeq(EntityPolySeq epolseq) {


	}

	@Override
	public void newExptl(Exptl exptl) {


	}

	@Override
	public void newCell(Cell cell) {

	}

	@Override
	public void newSymmetry(Symmetry symmetry) {

	}

	@Override
	public void newStructNcsOper(StructNcsOper sNcsOper) {

	}
	
	@Override
	public void newAtomSites(AtomSites atomSites) {
		
	}

	@Override
	public void newPdbxEntityNonPoly(PdbxEntityNonPoly pen) {


	}

	@Override
	public void newPdbxNonPolyScheme(PdbxNonPolyScheme ppss) {


	}

	@Override
	public void newPdbxPolySeqScheme(PdbxPolySeqScheme ppss) {


	}

	@Override
	public void newRefine(Refine r) {


	}

	@Override
	public void newStructAsym(StructAsym sasym) {


	}

	@Override
	public void newStructKeywords(StructKeywords kw) {


	}

	@Override
	public void newStructRef(StructRef sref) {


	}

	@Override
	public void newStructRefSeq(StructRefSeq sref) {


	}

	@Override
	public void newStructRefSeqDif(StructRefSeqDif sref) {


	}

	@Override
	public void setStruct(Struct struct) {


	}

	@Override
	public void newGenericData(String category, List<String> loopFields,
			List<String> lineData) {
		//System.out.println("unhandled category: " + category);

	}


	@Override
	public void newAuditAuthor(AuditAuthor aa)
	{


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


	}

	@Override
	public void newChemCompDescriptor(ChemCompDescriptor ccd) {
		ChemComp cc = dictionary.getChemComp(latestChemCompId);
		cc.getDescriptors().add(ccd);

	}

	@Override
	public void newPdbxStructOperList(PdbxStructOperList structOper) {


	}

	@Override
	public void newPdbxStrucAssembly(PdbxStructAssembly strucAssembly) {


	}

	@Override
	public void newPdbxStrucAssemblyGen(PdbxStructAssemblyGen strucAssembly) {


	}

	@Override
	public void newChemCompAtom(ChemCompAtom atom) {
		dictionary.getChemComp(latestChemCompId).getAtoms().add(atom);
	}

	@Override
	public void newPdbxChemCompIndentifier(PdbxChemCompIdentifier id) {


	}

	@Override
	public void newChemCompBond(ChemCompBond bond) {
		dictionary.getChemComp(latestChemCompId).getBonds().add(bond);
	}

	@Override
	public void newPdbxChemCompDescriptor(PdbxChemCompDescriptor desc) {


	}

	@Override
	public void newEntitySrcGen(EntitySrcGen entitySrcGen) {


	}
	@Override
	public void newEntitySrcNat(EntitySrcNat entitySrcNat) {


	}

	@Override
	public void newEntitySrcSyn(EntitySrcSyn entitySrcSyn) {


	}

	@Override
	public void newStructConn(StructConn structConn) {


	}

	@Override
	public void newStructSiteGen(StructSiteGen gen) {

	}

	@Override
	public void newStructSite(StructSite site) {

	}

	@Override
	public void newEntityPoly(EntityPoly entityPoly) {

		
	}

	@Override
	public void newPdbxAuditRevisionHistory(PdbxAuditRevisionHistory history) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void newPdbxDatabaseStatus(PdbxDatabaseStatus status) {
		// TODO Auto-generated method stub
		
	}
}

